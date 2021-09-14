# -*- coding: utf-8 -*-


from pathlib import Path
from molvs import validate_smiles
from molvs import Standardizer
from rdkit import Chem
import streamlit as st
import base64

#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Tools - LISTo',
    layout='wide')

######
# Function to put a picture as header   
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

from PIL import Image
image = Image.open('cropped-header.png')
st.image(image)


st.markdown("![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)")
st.subheader(":pushpin:" "About Us")
st.markdown("We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.")
st.markdown(":computer:""**Web Site** " "<https://lideb.biol.unlp.edu.ar>")

#---------------------------------#
st.write("""
# LIDeB Tools - LISTo

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

""")

text = '''
---

'''



# OPTIONS

#---------------------------------#
# Sidebar - Collects user input features into dataframe
st.sidebar.header('Upload your smiles')

uploaded_file_1 = st.sidebar.file_uploader("Upload your smiles in a TXT file", type=["txt"])

st.sidebar.markdown("""
[Example TXT input file](https://raw.githubusercontent.com/cbellera/LISTo/main/diclofenac_example.txt)
""")



LISTo_setting = st.sidebar.checkbox('Check to change the default configuration', value=False)

if LISTo_setting == True:
    
    st.sidebar.write("-----")
    remove_stereo = st.sidebar.checkbox("remove stereochemistry", value=True)
    st.sidebar.markdown("*All stereochemistry information is removed*")
    st.sidebar.write("-----")
    
    tautomerize = st.sidebar.checkbox("more stable tautomer", value=True)  
    st.sidebar.markdown("*Return the tautomer parent of a given molecule*")
    st.sidebar.write("-----")
    
    remove_charges = st.sidebar.checkbox("remove charges", value=True)
    st.sidebar.markdown("*Return the uncharged version of the fragment parent*")
    st.sidebar.write("-----")
    
    replace_isotopes = st.sidebar.checkbox("replace isotopes", value=True)
    st.sidebar.markdown("*All atoms are replaced with the most abundant isotope for that element*")
    st.sidebar.write("-----")
    
    add_hydrogens = st.sidebar.checkbox("add explicit hydrogens", value=False)
    st.sidebar.markdown("*Return the molecule with explicit hydrogens*")
    st.sidebar.write("-----")
else:
    pass

st.sidebar.title(":speech_balloon: Contact Us")
st.sidebar.info(
"""
If you are looking to contact us, please
[:e-mail:](mailto:lideb@biol.unlp.edu.ar) or [Twitter](https://twitter.com/LIDeB_UNLP)
""")

####---------------------------------------------------------------------------####
#### Standarization by MOLVS ####

import pandas as pd
 
def LISTo(uploaded_file_1):   
    data = pd.read_csv(uploaded_file_1,sep="\t",header=None,index_col=None)
    molecules = list(data[0])
    standarized_smiles_list = []
    log_validation = []
    t = st.empty()

    for i, line in enumerate(molecules,start = 1):
        smiles = line.strip()
        mol = Chem.MolFromSmiles(smiles)
        s = Standardizer()
        # To log the problems in SMILES
        result = validate_smiles(smiles)
        if result == []:
            log_validation.append(f"SMILES {i}: standardized!" )
        else:
            log_validation.append(f"SMILES {i}: {result}" )

        try:
            if LISTo_setting == False:
                # The standardization process consists of the following stages: RDKit RemoveHs(), SanitizeMol(), MetalDisconnector, Normalizer, Reionizer
                standarized_mol = s.fragment_parent(mol) # Standardized the molecule and return the fragment parent of a given molecule, the largest organic covalent unit in the molecule
                standarized_mol = s.stereo_parent(standarized_mol, skip_standardize= True) #Return The stereo parent of a given molecule, has all stereochemistry information removed from tetrahedral centers and double bonds.
                standarized_mol = s.tautomer_parent(standarized_mol, skip_standardize=True) # Return the tautomer parent of a given molecule.
                standarized_mol = s.charge_parent(standarized_mol, skip_standardize= True) #Return the charge parent of a given molecule,  the uncharged version of the fragment parent
                standarized_mol = s.isotope_parent(standarized_mol, skip_standardize= True) #Return the isotope parent of a given molecule, has all atoms replaced with the most abundant isotope for that element.    
                standarized_mol = Chem.rdmolops.RemoveHs(standarized_mol)
            else:
                # The standardization process consists of the following stages: RDKit RemoveHs(), SanitizeMol(), MetalDisconnector, Normalizer, Reionizer
                standarized_mol = s.fragment_parent(mol) # Standardized the molecule and return the fragment parent of a given molecule, the largest organic covalent unit in the molecule
                if remove_stereo:
                    standarized_mol = s.stereo_parent(standarized_mol, skip_standardize= True) #Return The stereo parentof a given molecule, has all stereochemistry information removed from tetrahedral centers and double bonds.
                if tautomerize:
                    standarized_mol = s.tautomer_parent(standarized_mol, skip_standardize=True) # Return the tautomer parent of a given molecule.
                if remove_charges:
                    standarized_mol = s.charge_parent(standarized_mol, skip_standardize= True) #Return the charge parent of a given molecule,  the uncharged version of the fragment parent
                if replace_isotopes:
                    standarized_mol = s.isotope_parent(standarized_mol, skip_standardize= True) #Return the isotope parent of a given molecule, has all atoms replaced with the most abundant isotope for that element.    
                if add_hydrogens:
                    standarized_mol = Chem.rdmolops.AddHs(standarized_mol)
                else:
                    standarized_mol = Chem.rdmolops.RemoveHs(standarized_mol)
                    
            smile_standar = Chem.MolToSmiles(standarized_mol)
            standarized_smiles_list.append(smile_standar)
        except:
            st.write(f'Molecule {i} could not be standardized')
        t.markdown("Progress: Molecule " + str(i) +"/" + str(len(molecules)))   

    standarized_smiles_final = pd.DataFrame(standarized_smiles_list)
    log_validation_final = pd.DataFrame(log_validation)
    return standarized_smiles_final, log_validation_final

# To export files

def filedownload(df):
    txt = df.to_csv(sep="\t",index=False,header=False)
    b64 = base64.b64encode(txt.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/txt;base64,{b64}" download="standardized_SMILES.txt">Download TXT File with your standardized SMILES</a>'
    return href

def filedownload1(df):
    txt = df.to_csv(sep="\t",index=False,header=False)
    b64 = base64.b64encode(txt.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/txt;base64,{b64}" download="log_SMILES.txt">Download TXT File with the log</a>'
    return href


if uploaded_file_1 is not None:
    run = st.button("Standardize")
    if run == True:
        st.markdown("**Standardizing...**")
        standarized_smiles_final, log_validation_final = LISTo(uploaded_file_1)
        
        st.markdown(":point_down: **Here you can dowload the standarized SMILES**", unsafe_allow_html=True)
        st.markdown(filedownload(standarized_smiles_final), unsafe_allow_html=True)
        
        st.markdown(":point_down: **Here you can dowload the log**", unsafe_allow_html=True)
        st.markdown(filedownload1(log_validation_final), unsafe_allow_html=True)
else:
    st.markdown("""
         ** :point_left: Please upload your smiles on the left **
         """)

    st.info('Awaiting for TXT file to be uploaded.')


#Footer edit

footer="""<style>
a:link , a:visited{
color: blue;
background-color: transparent;
text-decoration: underline;
}
a:hover,  a:active {
color: red;
background-color: transparent;
text-decoration: underline;
}
.footer {
position: fixed;
left: 0;
bottom: 0;
width: 100%;
background-color: white;
color: black;
text-align: center;
}
</style>
<div class="footer">
<p>Made in  🐍 and <img style='display: ; ' href="https://streamlit.io" src="https://i.imgur.com/iIOA6kU.png" target="_blank"></img> Developed with ❤️ by <a style='display: ; text-align: center' href="https://twitter.com/capigol" target="_blank">Lucas </a> , <a style='display: ; text-align: center' href="https://twitter.com/denis_prada" target="_blank">Denis </a> and <a style='display: ; text-align: center' href="https://twitter.com/carobellera" target="_blank">Caro </a> for <a style='display:; text-align: center;' href="https://lideb.biol.unlp.edu.ar/" target="_blank">LIDeB</a></p>
</div>
"""
st.markdown(footer,unsafe_allow_html=True)
