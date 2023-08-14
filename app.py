from Bio.Blast import NCBIWWW, NCBIXML
import time
import warnings
import streamlit as st
from PIL import Image

img = Image.open("image.png")

st.image(img, use_column_width=True)

st.write(""" ### Your results will display here :
""")
def set_button_style():
    button_style = """
        <style>
            .stButton button {
                background-color: #FF6F61;
                color: white;
                border-radius: 100px;
                border: 2px solid #FFFFFF; 
                box-shadow: 2px 2px 5px #000000;
                width: 160px;
            }
        </style>
    """
    st.markdown(button_style, unsafe_allow_html=True)

    
st.sidebar.title('Protein alignment with BLAST')
seq = st.sidebar.text_area('Enter a sequence to analyze : ')
type = st.sidebar.text_input('Enter BLAST type : ',placeholder='exemple : blastp, blastn ...')
data = st.sidebar.text_input('Enter the database : ', placeholder='exemple : nt, nr ...')
detailed = st.sidebar.text_input('Enter detailed search : ', placeholder='Homo sapiens')
set_button_style()
b = st.sidebar.button("Run BLAST")
s = st.sidebar.button("Run detailed search")



def blast(seq, type, data):

    st.write("""Starting BLAST search...""") 
    start_time = time.time()
    result_handle2 = NCBIWWW.qblast(type, data, seq)
    end_time = time.time()
    
    st.write("""BLAST search completed.""")
    st.write(f"Time taken: {end_time - start_time:.2f} seconds.")
    
    with open("output.xml", "w") as out_handle:
        out_handle.write(result_handle2.read())
    st.write("Results saved to output.xml")
    st.write("===================================")
    result_handle = open("output.xml")
    blast_record = NCBIXML.read(result_handle)
    st.write(""" **Results :** """)
    count = 0
    for alignment in blast_record.alignments:
        st.write("**"+str(alignment.title)+"**")
        for hsp in alignment.hsps:
            st.write(str(hsp.query[0:50]))
            st.write(str(hsp.query[0:50]))
            st.write(str(hsp.query[0:50]))
            percentage_result = (hsp.query_end - hsp.query_start + 1) / hsp.align_length * 100
            st.write(str(hsp))
            if percentage_result >= 50:  
                st.write(f'<p style="color:green;">{percentage_result:.2f}</p>', unsafe_allow_html=True)
            else:
                st.write(f'<p style="color:red;">{percentage_result:.2f}</p>', unsafe_allow_html=True)
        count += 1
        if count == 10:
            break

def detailed_search(detailed):
    
    result_handle = open("output.xml")
    blast_record = NCBIXML.read(result_handle)
    
    st.write(""" Detailed search on : """ + "**" + str(detailed) + "**")
    for alignment in blast_record.alignments:
        if detailed in alignment.title:
            st.write(str(alignment.title))
            for hsp in alignment.hsps:
                st.write(str(hsp.query[0:50]))
                st.write(str(hsp.match[0:50]))
                st.write(str(hsp.sbjct[0:50]))
                st.write(str(hsp.positives) + str(hsp.score) + str(hsp.expect) + str(hsp.bits) + str(hsp.gaps))

def save_seq(seq):
    with open("sequenceinput.fsa", "w") as file:
        file.write(str(seq))
    st.write("sequence saved to sequenceinput.fsa")

if seq :
    save_seq(seq)

if seq and  type  and data and b :
    blast(seq, type, data)

if detailed and s : 
    detailed_search(detailed)


