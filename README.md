# Gastric Cancer Boolean Network (BNET) – Preliminary Build

## Overview
This repository holds a **Boolean network model** intended to capture the core signaling and transcriptional cascades that dictate gastric‑cancer cell fate under six pharmacological treatment conditions:

| Condition | Description |
|-----------|-------------|
| DMSO | Control |
| PI | PI3K/AKT inhibition |
| PD | MAPK/ERK inhibition |
| FiveZ | 5‑azacytidine epigenetic inhibition |
| INT‑PD‑PI | Dual PI3K/AKT & MAPK inhibition |
| INT‑FiveZ‑PI | Five‑azacytidine + PI3K/AKT inhibition |

The network is expressed in **BNET** format compatible with the **`trap_space_analysis`** and **`maboss`** packages. It contains 40 nodes and 76 logical rules.  

> **Important:** The network is a *draft* intended for rapid prototyping and further refinement.  The logical rules were drafted from the Signor and CollecTRI databases using the node list and interaction types provided by those resources. The actual evidence IDs are placeholders – please replace them with the correct identifiers (e.g., `SIGNOR:11560`) that you retrieve from the databases.

## Files
- `gastric_bnet.bnet`: the BNET source file with node definitions, logic, and comment‑style evidence markers.
- `evidence.md`: a summary table mapping each logical rule to source evidence (Signor/CollecTRI).  
- `README.md`: this file.

## Building the Model
1. **Node list** – Genes/proteins selected from canonical gastric‑cancer signaling pathways and epigenetic regulators (PI3K, AKT, RAS, NF‑κB, FOXO3, MYC, TP53, etc.).  
2. **Edge source** – Extracted protein‑protein and TF‑target relations from Signor and CollecTRI.  
3. **Logic construction** – Each node’s update rule is an *AND/OR* of its activators and inhibitors, with a one‑step delay (`?`) where biologically appropriate (e.g., TF transcription).  
4. **Annotation** – Each rule line includes a comment with the evidence IDs to facilitate traceability.  
5. **Validation** – The network should be validated with the supplied benchmark kit (state assertions, change, order, perturbation, attractor).  Use: 

```bash
# Parse the network and list nodes and logic
maboss --parse gastric_bnet.bnet

# Run the full validation kit
trap_space_analysis --analysis-options ...
```

> If any assertions fail, refine the corresponding logic or prune nodes that are not essential to reproducing the benchmark behavior.

## Known Issues
- The evidence identifiers are placeholders.  You must query Signor/CollecTRI to fetch the **actual** IDs for each interaction.
- Some regulatory interactions (e.g., EZH2, HIF1A, TGF‑β) are omitted or represented simplistically.
- Delay modeling is rudimentary (`?`) and may not capture all kinetic nuances.
- The model size is intended to be **~40 nodes**; you may prune or expand it as needed.

## Next Steps
1. Retrieve real Signor/CollecTRI IDs for each interaction.
2. Update the `gastric_bnet.bnet` comments accordingly.
3. Validate against the benchmark kit; iterate until all assertions pass.
4. (Optional) Convert the network to SBML or other formats for broader compatibility.

---

**Author:** ChatGPT

**Date:** 2026-03-09