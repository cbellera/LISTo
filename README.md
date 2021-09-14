📌About Us

We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.

💻Web Site https://lideb.biol.unlp.edu.ar


## LIDeB Tools - LISTo

LIDeB's Standardization Tool is a WebApp to standardize SMILES based in MOLVS.

The tool uses the following packages: [RDKIT](https://www.rdkit.org/docs/index.html),  [MOLVS](https://molvs.readthedocs.io/)

Default setting will perform next actions to each smiles:

- remove explicit hydrogens by RemoveHs (RDKIT)
- Kekulize, check valencies, set aromaticity, conjugation and hybridization by SanitizeMol (RDKIT)
- disconnect metal atoms that are defined as covalently bonded to non-metals by MetalDisconnector (RDKIT)
- correct functional groups and recombining charges by Normalizer (RDKIT)
- fix charges and reionize a molecule by Reionizer (RDKIT)
- select the largest covalent fragment by fragment_parent (molvs)
- remove stereochemistry by stereo_parent (molvs)
- select the more stable tautomeric form by tautomer_parent (molvs)
- remove charges - neutralization - by charge_parent (molvs)
- select the most abundant isotope by isotope_parent (molvs)


If you are looking to contact us, please send a mail to lideb@biol.unlp.edu.ar or contact us by Twitter (https://twitter.com/LIDeB_UNLP)
