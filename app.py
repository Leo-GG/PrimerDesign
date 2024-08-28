import streamlit as st
from PrimerUtils import utils
import requests
import pandas as pd

# Function to run with user input
def process_inputs(chromosome, target_position, region_length=100, primer_length=20, target_tm=60, min_gc=40, max_gc=60):
    
    start= target_position - region_length // 2 - primer_length
    end = target_position + region_length // 2 + primer_length
    
    try:
        ensembl_response = utils.get_genome_sequence(chromosome, start, end)
        sequence=ensembl_response['seq']
        build=ensembl_response['id'].split(':')[1]

        result={ 
            'Build':build,
            'Sequence':sequence
        }
        
    except requests.exceptions.HTTPError as e:
        # Handle the HTTPError separately
        #print(f"HTTPError: {e}")
        result={'HTTPError': e}
        return result
        
    
    except requests.exceptions.RequestException as e:
        # Handle other possible exceptions related to the request (e.g., ConnectionError)
        #print(f"RequestException: {e}")
        result={'RequestException': e}
        return result
    
    if ((len(sequence)<(2*primer_length+region_length)) | ('N' in sequence)) :
        result['Error']=f"Target sequence at CHRM {chromosome} POS {target_position} with given length is not suitable"
        return result
    
    try:
        
        forward_primer, reverse_primer, forward_tm, reverse_tm = utils.design_primers(sequence,
                                                                                      primer_length=primer_length,
                                                                                      target_tm=target_tm, 
                                                                                      min_gc=min_gc, 
                                                                                      max_gc=max_gc)

        result['Forward primer']=forward_primer
        result['Reverse primer']=reverse_primer
        result['Forward primer Tm']=f"{forward_tm:.2f}°C"
        result['Reverse primer Tm']=f"{reverse_tm:.2f}°C"
        
        try:
            result['Blast results']=utils.run_in_silico_pcr(forward_primer,reverse_primer)
            
        except requests.exceptions.HTTPError as e:
            # Handle the HTTPError separately
            #print(f"HTTPError: {e}")
            result['Blast HTTPError']=e
            return result   
        except requests.exceptions.RequestException as e:
            # Handle other possible exceptions related to the request (e.g., ConnectionError)
            #print(f"RequestException: {e}")
            result['Blast RequestException']=e
            return result
        
        return result

    except ValueError as e:
        
        result['Error'] = str(e)
        
        return result

# Streamlit app
def main():
    st.title("Quick primer design")
    description = r""" This application creates PCR primers to amplify the sequence around a selected position in the **human** genome. It returns the 
    pair of forward and reversed primers and then uses the [UCSC In-Silico PCR tool](https://genome.ucsc.edu/cgi-bin/hgPcr) to check for
    alignments along the genome.
    
    Tm are calculated using the [MeltingTemp](https://biopython.org/docs/1.75/api/Bio.SeqUtils.MeltingTemp.html) module from the [Biopython](https://biopython.org/) package.
    """
    st.markdown(description, unsafe_allow_html=True)


    # Sidebar (left panel)
    st.sidebar.title("Input Panel")

    # Dropdown for choosing chromosome
    chromosome_list=list(['CHRM '+str(i) for i in range(1,23)])+list(['CHRM X','CHRM Y'])
    chromosome = st.sidebar.selectbox(
        "Choose Chromosome",
        chromosome_list,0
    ).split('CHRM ')[1]
        
    
    # Input box for position (with error handling for numerical input)
    position = st.sidebar.text_input("Enter Position",1)
    try:
        position = int(position)
    except ValueError:
        st.sidebar.error("Position must be a number")
        position = None
    
    # Radio buttons to choose organism
    organism = st.sidebar.radio(
        "Choose Organism",
        (["Human"])
    )
    
    # Additional numerical inputs for extra parameters
    region_length = st.sidebar.number_input("Target region length", value=100)
    primer_length = st.sidebar.number_input("Initial primer length", value=20)
    target_tm = st.sidebar.number_input("Ideal Tm", value=55)
    min_gc = st.sidebar.number_input("Minimum %GC", value=20)
    max_gc = st.sidebar.number_input("Maximum %GC", value=80)
    
    # Run button
    run_button = st.sidebar.button("Run")
    
    # Main page
    

    if run_button:
        st.title("Results")
        # Only run if position is a valid number
        result = process_inputs(chromosome, position, region_length, primer_length, target_tm, min_gc, max_gc)
        if position is not None:
            try:
                # Call the processing function with the selected parameters
                result = process_inputs(chromosome, position, region_length, primer_length, target_tm, min_gc, max_gc)
                # Display the results
                st.markdown("### Primers")
                for key, value in result.items():
                    if ('Blast' in key):
                        blast_output=value
                    else:
                        st.markdown(f"**{key}:** {value}")
            except:
                st.write(result)
                
            st.markdown("### Blast search")
            df=pd.DataFrame(blast_output)
            st.dataframe(df)
            
        else:
            st.write("Please enter a valid position.")
            
        
        
    else:
        st.title("Press Run to design your primers")
    
if __name__ == "__main__":
    main()