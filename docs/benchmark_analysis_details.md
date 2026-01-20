# Benchmark Analysis Details: Tool Equivalence and Metrics

This document details the configuration equivalence between **Seqproc**, **Matchbox**, and **Splitcode** for the paper benchmarks, and clarifies the specific metrics used for "validity" and "recovery."

## 1. Tool Equivalence in Chemistry Parsing

All three tools are configured to parse the same biological structures, but they use different abstractions (Geometry, Pattern Matching, Tag Graphs) to achieve this.

### A. Sci-Seq 3 Chemistry
**Target Structure (Read 1):** `[BC1: 9-10bp] - [Anchor: CAGAGC] - [UMI: 8bp] - [BC2: 10bp]`

The primary challenge is the **variable length** of BC1 (9 or 10bp), which shifts the position of the Anchor and subsequent elements.

| Feature | **Seqproc** (`sciseq3.geom`) | **Matchbox** (`sciseq3.mb`) | **Splitcode** (`sciseq3.config`) |
| :--- | :--- | :--- | :--- |
| **Variable BC1** | **`brc1 = norm(b[9-10])`**<br>Defines a normalized barcode region that acts as a variable-length window (9 or 10 bases). | **Explicit Branching**<br>Defines two separate match patterns:<br>1. `[ bc1:|9| ... ]`<br>2. `[ bc1:|10| ... ]`<br>The tool attempts to match both. | **`@extract <...[9-10]>{{anchor}}`**<br>Extracts the sequence *preceding* the anchor. The `[9-10]` notation tells Splitcode the expected length range of that preceding segment. |
| **Anchor** | **`hamming(f[CAGAGC], 1)`**<br>Searches for `CAGAGC` with a maximum Hamming distance of **1**. | **`primer~0.17`**<br>Matches `CAGAGC` with ~17% error tolerance.<br>$(6 \text{bp} \times 0.17 \approx 1.02)$<br>Effective Hamming distance: **1**. | **`anchor CAGAGC 1`**<br>Defines the tag `CAGAGC` with an explicit mismatch distance of **1**. |
| **UMI / BC2** | **Positional Definition**<br>`umi = u[8]`, `brc2 = b[10]`<br>Defined to occur immediately after the anchor. | **Positional Pattern**<br>`umi:|8| bc2:|10|`<br>Defined to occur immediately after the primer in the match pattern. | **Relative Extraction**<br>`@extract {{anchor}}8<...[10]>`<br>Extracts 8bp (UMI) immediately following the anchor, then the next 10bp (BC2). |

**Equivalence Verdict:**
All three tools are functionally equivalent for Sci-Seq 3. They all dynamically locate the `CAGAGC` anchor allowing 1 mismatch and adjust the BC1 extraction boundaries accordingly.

---

### B. SPLiT-seq Paired-End (Replacement)
**Target Structure (Read 2):** `[NN] - [UMI: 8bp] - [BC3: 8bp] - [Linker 1: 30bp] - [BC2: 8bp] - [Linker 2: 30bp] - [BC1: 8bp]`

The challenge here is **whitelist correction** (mapping noisy sequences to a known list) and navigating long linkers.

| Feature | **Seqproc** (`splitseq_replacement.geom`) | **Matchbox** (`splitseq_replacement.mb`) | **Splitcode** (`splitseq_paper.config`) |
| :--- | :--- | :--- | :--- |
| **Linker 1 (30bp)** | **`hamming(f[GTGG...], 2)`**<br>Strict Hamming distance $\le 2$. | **`l1~0.2`**<br>Error rate $\approx 20\%$.<br>$(30 \text{bp} \times 0.2 = 6)$<br>Effective Hamming distance $\le 6$. | **`linker1 GTGG... 2`**<br>Strict Hamming distance $\le 2$. |
| **Linker 2 (30bp)** | **`hamming(f[ATCC...], 2)`**<br>Strict Hamming distance $\le 2$. | **`l2~0.2`**<br>Matches **16bp prefix** (`ATCC...`) with $\approx 20\%$ error ($\approx 3$ mismatches). | **`linker2 ATCC... 2`**<br>Strict Hamming distance $\le 2$. |
| **Barcode Matching** | **`map_with_mismatch(..., 2)`**<br>Maps against whitelist with max **2** mismatches. | **`bc3.round_23~0.13`**<br>Fuzzy match against CSV with $\approx 13\%$ error.<br>$(8 \text{bp} \times 0.13 \approx 1.04)$<br>Effective mismatches: **1**. | **`bc_X ... 2`**<br>Explicitly defines every valid barcode with distance **2**. |

**Equivalence Verdict:**
There is a **configuration discrepancy** in the current benchmark setup for SPLiT-seq:
1.  **Linkers:** Matchbox is configured significantly **looser** (allows ~6 mismatches) than Seqproc/Splitcode (allows 2).
2.  **Barcodes:** Matchbox is configured **stricter** (allows ~1 mismatch) than Seqproc/Splitcode (allows 2).

This explains why Matchbox might behave differently on low-quality reads. It will find more linkers (looser) but reject more barcodes (stricter).

---

## 2. Q&A: Detailed Explanations

### Q1: Why is Matchbox configured to be "looser" (or different) for SPLiT-seq?
**Short Answer:** The current Matchbox configuration uses error rates (`~0.2`) rather than fixed Hamming distances, leading to differences in strictness.

**Detailed Explanation:**
Matchbox uses a "fuzzy search" algorithm based on error rates (Levenshtein/Edit distance approximation), whereas Seqproc and Splitcode typically operate on fixed Hamming distances (substitutions only).

In the current `splitseq_replacement.mb` config:
*   `l1~0.2`: 20% error on a 30bp sequence is $30 \times 0.2 = 6$ edits. This is **much looser** than the Hamming-2 used by Seqproc/Splitcode.
*   `bc3.round_23~0.13`: 13% error on 8bp is $8 \times 0.13 \approx 1.04$ edits. This effectively allows **1 mismatch**, which is **stricter** than the Hamming-2 allowed by Seqproc/Splitcode.

**Why this discrepancy?**
This likely stems from the default examples or recommended settings for Matchbox, which often favor relative error rates over absolute distances. To make them strictly equivalent to the paper's definition (Hamming 2), the Matchbox config would need to be tuned carefully (e.g., setting specific error rates that strictly calculate to 2 edits, though Matchbox's model is probabilistic/approximate in some modes).

### Q2: Are the "valid reads" a percentage of Input or Output?
**Short Answer:** The "Recovery Rate" reported in the JSON/Figures is **% of Total Input**. The "Valid" counts discussed in the logs are the **intersection** of the tool's output and the ground-truth valid list.

**The "Skewed Validity" Problem:**
You correctly identified a potential pitfall: *"If we have lower recovery, then the similar valid percentage is kinda skewed."*

Let's define the metrics strictly:

1.  **Recall (Recovery Rate):**
    $$ \frac{\text{Reads Output by Tool}}{\text{Total Input Reads}} $$
    *   This is the main number in your bar charts.
    *   Example: Seqproc extracts 74%, Matchbox extracts 83%.

2.  **Precision (Validity of Output):**
    $$ \frac{\text{Valid Reads Output by Tool}}{\text{Total Reads Output by Tool}} $$
    *   This tells you: "Of the reads the tool decided to keep, how many were actually correct?"
    *   If Tool A is very conservative (low recovery) but only outputs perfect matches, it has **High Precision, Low Recall**.
    *   If Tool B is very loose (high recovery) but outputs lots of junk, it has **Low Precision, High Recall**.

**In our Benchmark:**
*   **Seqproc (SPLiT-seq PE):** ~74% Recovery.
*   **Matchbox (SPLiT-seq PE):** ~83% Recovery.

The fact that Matchbox has **higher** recovery (83% vs 74%) AND is finding valid reads (based on the validity checker) suggests it is genuinely recovering more data, likely due to the **looser Linker 1 setting** (allowing 6 mismatches vs 2). The "stricter" barcode setting (1 mismatch) suggests that the reads it *is* losing (compared to if it had dist-2) are outweighed by the reads it is *gaining* from the looser linker search.

**Conclusion:**
Matchbox is recovering **more total reads** because its linker search is more permissive. If we tightened the linker search to `~0.06` (approx 2/30), its recovery would likely drop to be closer to Seqproc/Splitcode.
