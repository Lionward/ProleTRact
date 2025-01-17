import pandas as pd 
import numpy as np


def count_motifs_with_zeros(motif_ids, num_motifs):
    # Use np.unique to get the unique values and their counts
    unique, counts = np.unique(motif_ids, return_counts=True)

    # Initialize an array of zeros for the full range of motifs
    full_counts = np.zeros(num_motifs, dtype=int)

    # Place the counts of the unique motifs in the appropriate positions
    full_counts[unique] = counts

    return tuple(full_counts)

def parse_vcf_tandemtwister(vcf, add_motif_counts=False):
    df = []
    c = 0
    for record in vcf.fetch():
        c += 1
        ID = record.id
        chrom = ID.split(':')[0]
        start = int(ID.split(':')[1].split('-')[0])
        
        end = int(ID.split(':')[1].split('-')[1])
        motifs = record.info['MOTIFS']
        if (type(motifs) == str):
            # make it as tuple
            motifs = tuple(motifs.split(','))
        
        # motif_lens = [len(s) for s in motifs]

        motif_ids_H1 = record.info['MOTIF_IDs_H1']
        motif_ids_H2 = record.info['MOTIF_IDs_H2']

        CN_H1 = int(record.info['CN_H1'])
        CN_H2 = int(record.info['CN_H2'])
        CN_ref = int(record.info['CN_ref'])

        if motif_ids_H2 == None and CN_H2 != 0 :
            motif_ids_H2 = motif_ids_H1
        
        # print(f'motif_ids_H1 = {motif_ids_H1}, type = {type(motif_ids_H1)}, len = {len(motif_ids_H1)}')
        # print(f'motif_ids_H2 = {motif_ids_H2}, type = {type(motif_ids_H2)}, len = {len(motif_ids_H2)}')
        

        len_h1, len_h2, motif_concat_h1, motif_concat_h2 = -1, -1, "NA", "NA"

        if CN_H1 > 0:
            len_h1 = sum([len(motifs[int(m_id)]) for m_id in motif_ids_H1])
            motif_concat_h1 = ''.join([motifs[int(m_id)] for m_id in motif_ids_H1])

        if CN_H2 > 0:
                
            len_h2 = sum([len(motifs[int(m_id)]) for m_id in motif_ids_H2])
            motif_concat_h2 = ''.join([motifs[int(m_id)] for m_id in motif_ids_H2])

        if add_motif_counts:
            motif_count_h1, motif_count_h2 = np.nan, np.nan
        
        if add_motif_counts and motif_ids_H1 != None and motif_ids_H2 != None:
            # motif counts
            ids_h1 = list(map(int, motif_ids_H1))
            ids_h2 = list(map(int, motif_ids_H2))
            num_motifs = len(motifs) #max(ids_h1 + ids_h2) + 1
            motif_count_h1 = count_motifs_with_zeros(ids_h1, num_motifs)
            motif_count_h2 = count_motifs_with_zeros(ids_h2, num_motifs)


        # get the GT value from the record string 
        GT = str(record).split('\t')[9].split(':')[0]

        if add_motif_counts:
            df.append([chrom, start, end, motifs, motif_ids_H1, motif_ids_H2, CN_H1, CN_H2, CN_ref, GT, len_h1, len_h2, motif_count_h1, motif_count_h2, motif_concat_h1, motif_concat_h2])
        else:
            df.append([chrom, start, end, motifs, motif_ids_H1, motif_ids_H2, CN_H1, CN_H2, CN_ref, GT, len_h1, len_h2, motif_concat_h1, motif_concat_h2])

    # convert to dataframe
    df = pd.DataFrame(df, columns=['chrom', 'start', 'end', 'motifs', 'motif_ids_H1', 'motif_ids_H2', 'CN_H1', 'CN_H2', 'CN_ref', 'GT', 'len_h1', 'len_h2', 'motif_concat_h1', 'motif_concat_h2'])
    return df

def parse_vcf_TRGT(vcf, add_motif_counts=False):
    df = []

    for record in vcf.fetch():
        ID = record.info['TRID']
        chrom = ID.split('_')[0]
        start = int(ID.split('_')[1])
        end = int(ID.split('_')[2])
        motifs = record.info['MOTIFS']        
        # make a list of the motifs from the tuple since tuple has no split method
        for sample in record.samples.values():
            mc_value = sample.get('MC')
            if (mc_value[0] == ".")  :
                continue

            # make list of the MC values[0]
            mc_values_hap1 = mc_value[0]
            if len(mc_value) == 1:
                mc_values_hap2 = mc_values_hap1
            else:
                if mc_value[1] == ".":
                    mc_values_hap2 = mc_values_hap1
                else:
                    mc_values_hap2 = mc_value[1]
                    
            mc_values_hap1 = list(map(int, mc_values_hap1.split('_')))
            mc_values_hap2 = list(map(int, mc_values_hap2.split('_')))
         
            # sort the motifs but sort the MC values of mc_values_hap1, and mc_values_hap2 with the same order, but make the sorting based on the motifs
            mc_values_hap1 = [x for _, x in sorted(zip(motifs, mc_values_hap1), key=lambda pair: pair[0])]
            mc_values_hap2 = [x for _, x in sorted(zip(motifs, mc_values_hap2), key=lambda pair: pair[0])]
  
            
            # update the motifs with the sorted motifs
            motifs = sorted(motifs)
            # make the motfs as tuple 
            motifs = tuple(motifs)


            motifs_ids_hap1 = [str(mc_values_hap1.index(i))*i for i in mc_values_hap1]
            motifs_ids_hap2 = [str(mc_values_hap2.index(i))*i for i in mc_values_hap2]
            motifs_ids_hap1 = [char for string in motifs_ids_hap1 for char in string]
            motifs_ids_hap2 = [char for string in motifs_ids_hap2 for char in string]
            # make them tuple again
            motifs_ids_hap1 = tuple(motifs_ids_hap1)
            motifs_ids_hap2 = tuple(motifs_ids_hap2)
            GT = sample['GT']
   
            # make the tuple like 1/2 and so on
            GT = '/'.join([str(i) for i in GT])
            CN_hap1 = len(motifs_ids_hap1)
            CN_hap2 = len(motifs_ids_hap2)
            CN_ref = np.nan
            if GT == '0/0' or GT == '1/1':
                CN_hap2 = CN_hap1
                motifs_ids_hap2 = motifs_ids_hap1
            

            # print(f'motif_ids_H1 = {motifs_ids_hap1}, type = {type(motifs_ids_hap1)}, len = {len(motifs_ids_hap1)}')
            # print(f'motif_ids_H2 = {motifs_ids_hap2}, type = {type(motifs_ids_hap2)}, len = {len(motifs_ids_hap2)}')

            len_h1, len_h2, motif_concat_h1, motif_concat_h2 = -1, -1, "NA", "NA"

            if CN_hap1 > 0:
                len_h1 = sum([len(motifs[int(m_id)]) for m_id in motifs_ids_hap1])
                motif_concat_h1 = ''.join([motifs[int(m_id)] for m_id in motifs_ids_hap1])

            if CN_hap2 > 0:
                len_h2 = sum([len(motifs[int(m_id)]) for m_id in motifs_ids_hap2])
                motif_concat_h2 = ''.join([motifs[int(m_id)] for m_id in motifs_ids_hap2])
            
            if add_motif_counts:
                motif_count_h1, motif_count_h2 = np.nan, np.nan
            
            if add_motif_counts and motifs_ids_hap1 != None and motifs_ids_hap2 != None:
                # # motif counts
                # ids_h1 = list(map(int, motifs_ids_hap1))
                # ids_h2 = list(map(int, motifs_ids_hap2))
                # # print('type of ids =', type(ids_h1), type(ids_h2))
                # num_motifs = len(motifs) # max(ids_h1 + ids_h2) + 1
                # motif_count_h1 = count_motifs_with_zeros(ids_h1, num_motifs)
                # motif_count_h2 = count_motifs_with_zeros(ids_h2, num_motifs)


                ids_h1 = np.array(motifs_ids_hap1, dtype=int)
                ids_h2 = np.array(motifs_ids_hap2, dtype=int)
                
                # Check and print the type if it's not integer
                if ids_h1.dtype != int:
                    print('ids_h1:', ids_h1, 'Type:', ids_h1.dtype)
                    
                if ids_h2.dtype != int:
                    print('ids_h2:', ids_h2, 'Type:', ids_h2.dtype)
                
                num_motifs = len(motifs)  # max(ids_h1 + ids_h2) + 1
                motif_count_h1 = count_motifs_with_zeros(ids_h1, num_motifs)
                motif_count_h2 = count_motifs_with_zeros(ids_h2, num_motifs)
                
            if add_motif_counts: 
                df.append([chrom, start, end, motifs, motifs_ids_hap1, motifs_ids_hap2, CN_hap1, CN_hap2, CN_ref, GT, len_h1, len_h2, motif_count_h1, motif_count_h2, motif_concat_h1, motif_concat_h2])
            else:
                df.append([chrom, start, end, motifs, motifs_ids_hap1, motifs_ids_hap2, CN_hap1, CN_hap2, CN_ref, GT, len_h1, len_h2, motif_concat_h1, motif_concat_h2])
    
    # convert to dataframe
    df = pd.DataFrame(df, columns=['chrom', 'start', 'end', 'motifs', 'motif_ids_H1', 'motif_ids_H2', 'CN_H1', 'CN_H2', 'CN_ref', 'GT', 'len_h1', 'len_h2', 'motif_concat_h1', 'motif_concat_h2'])
    
    return df


def parse_vcf_vamos(vcf, add_motif_counts=False):
    df = []
    for record in vcf.fetch():
        chrom = record.chrom
        start = record.start +1
        end = record.stop
        motifs = record.info['RU']
        
        # check if ALTANNO_H1 is in the record
        if 'ALTANNO_H1' in record.info:
            motif_ids_H1 = record.info['ALTANNO_H1']
            
        else:
            motif_ids_H1 = ""
        
        if 'ALTANNO_H2' in record.info:
            motif_ids_H2 = record.info['ALTANNO_H2']
           
        else:
            motif_ids_H2 = ""

        CN_H1 = int(record.info['LEN_H1'])
        if 'LEN_H2' in record.info:
            CN_H2 = int(record.info['LEN_H2'])
        else:
            CN_H2 = CN_H1
            motif_ids_H2 = motif_ids_H1

        CN_ref = np.nan
        GT = record.samples[0]['GT']

        # print(f'motif_ids_H1 = {motif_ids_H1}, type = {type(motif_ids_H1)}, len = {len(motif_ids_H1)}')
        # print(f'motif_ids_H2 = {motif_ids_H2}, type = {type(motif_ids_H2)}, len = {len(motif_ids_H2)}')

        len_h1, len_h2, motif_concat_h1, motif_concat_h2 = -1, -1, "NA", "NA"

        if CN_H1 > 0:
            len_h1 = sum([len(motifs[int(m_id)]) for m_id in motif_ids_H1])
            motif_concat_h1 = ''.join([motifs[int(m_id)] for m_id in motif_ids_H1])

        if CN_H2 > 0:
            len_h2 = sum([len(motifs[int(m_id)]) for m_id in motif_ids_H2])
            motif_concat_h2 = ''.join([motifs[int(m_id)] for m_id in motif_ids_H2])
        
        if add_motif_counts:
            motif_count_h1, motif_count_h2 = np.nan, np.nan
        
        if add_motif_counts and motif_ids_H1 != None and motif_ids_H2 != None:
            # motif counts
            ids_h1 = list(map(int, motif_ids_H1))
            ids_h2 = list(map(int, motif_ids_H2))
            num_motifs = len(motifs) # max(ids_h1 + ids_h2) + 1
            motif_count_h1 = count_motifs_with_zeros(ids_h1, num_motifs)
            motif_count_h2 = count_motifs_with_zeros(ids_h2, num_motifs)


        GT = '/'.join([str(i) for i in GT])


        if add_motif_counts:
            df.append([chrom, start, end, motifs, motif_ids_H1, motif_ids_H2, CN_H1, CN_H2, CN_ref, GT, len_h1, len_h2, motif_count_h1, motif_count_h2, motif_concat_h1, motif_concat_h2])
        else:
            df.append([chrom, start, end, motifs, motif_ids_H1, motif_ids_H2, CN_H1, CN_H2, CN_ref, GT, len_h1, len_h2, motif_concat_h1, motif_concat_h2])


    # convert to dataframe
    df = pd.DataFrame(df, columns=['chrom', 'start', 'end', 'motifs', 'motif_ids_H1', 'motif_ids_H2', 'CN_H1', 'CN_H2', 'CN_ref', 'GT', 'len_h1', 'len_h2', 'motif_concat_h1', 'motif_concat_h2'])
    return df