import requests
import sys
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from bs4 import BeautifulSoup


def get_genome_sequence(chromosome, start, end):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1"
    headers = {"Content-Type": "application/json"}
 
    
    response = requests.get(server + ext, headers=headers)
    # Raise an HTTPError for bad responses (4xx and 5xx)
    response.raise_for_status()


    #if not response.ok:
    #    response.raise_for_status()
    #    sys.exit()
    return response.json()    
    


def design_primers(sequence, primer_length=20, target_tm=60, min_gc=40, max_gc=60, max_iterations=10):
    forward_primer = sequence[0:primer_length]
    reverse_primer =sequence[len(sequence)-primer_length:len(sequence)]
 
    # Get reverse complement for the reverse primer
    reverse_complement = str(Seq(reverse_primer).reverse_complement())
    
    # Adjust primer length to achieve target Tm (with a maximum number of iterations)
    for _ in range(max_iterations):
        f_tm = calculate_tm(forward_primer)
        r_tm = calculate_tm(reverse_complement)
        
        if abs(f_tm - target_tm) <= 2 and abs(r_tm - target_tm) <= 2:
            break
        
        if f_tm < target_tm and len(forward_primer) < 30:
            forward_primer += sequence[len(forward_primer)]
        elif f_tm > target_tm and len(forward_primer) > 18:
            forward_primer = forward_primer[:-1]
       
        if r_tm < target_tm and len(reverse_complement) < 30:
            reverse_complement =  reverse_complement + str(Seq(sequence[-len(reverse_complement)-1]).complement())
        elif r_tm > target_tm and len(reverse_complement) > 18:
            reverse_complement = reverse_complement[1:]
    
    # Check GC content (40-60% is ideal)
    forward_gc = (forward_primer.count('G') + forward_primer.count('C')) / len(forward_primer) * 100
    reverse_gc = (reverse_complement.count('G') + reverse_complement.count('C')) / len(reverse_complement) * 100
    
    if min_gc <= forward_gc <= max_gc and min_gc <= reverse_gc <= max_gc:
        return forward_primer, reverse_complement, f_tm, r_tm
    else:
        raise ValueError(f"GC content is out of the ideal range ({min_gc}-{max_gc}%)")
        
def calculate_tm(primer):
    return mt.Tm_NN(primer, nn_table=mt.DNA_NN4)


def run_in_silico_pcr(forward_primer, reverse_primer, genome_assembly='hg38'):
    url = 'https://genome.ucsc.edu/cgi-bin/hgPcr'
    payload = {
        'org': 'Human',
        'db': 'hg38',
        'wp_target': 'genome',
        'wp_f': forward_primer,
        'wp_r': reverse_primer,
        'wp_size': '4000',
        'Submit': 'Submit',
        'wp_perfect':'5',
        'wp_good':'15',
        'boolshad.wp_flipReverse':'0',
        'wp_append':'on',
        'boolshad.wp_append':'0'             
    }
    
    response = requests.post(url, data=payload)
    
    response.raise_for_status()

    if response.status_code == 200:
        return extract_alignments(response.text)
    else:
        return [{'Blast error':response.status_code}]
    

def extract_alignments(html_text):
    soup = BeautifulSoup(html_text, 'html.parser')
    alignments = []

    for pre in soup.find_all('pre'):
        results_= pre.text.strip().split('>')[1:]     
        for r in results_:
            lines = r.split('\n')
            alignment={}
            for line in lines:
                if line.startswith('chr'):
                    parts = line.split()
                    location = parts[0]
                    length = parts[1].replace('bp', '')
                    forward_primer = parts[2]
                    reverse_primer = parts[3]
                    sequence = '' #'\n'.join(parts[4:]) if len(parts) > 4 else ''
                    alignment={'location': location,
                        'length': length,
                        'forward_primer': forward_primer,
                        'reverse_primer': reverse_primer,
                        'sequence': sequence}
                elif len(alignment)>0:
                    alignment['sequence']='\n'.join([alignment['sequence'],line])
                    
                alignments.append(alignment)
    
    return alignments