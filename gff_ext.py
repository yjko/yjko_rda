'''
Created on 2016. 7. 13.

@author: yjko
'''

fsl = 100000

sp_list = ['arabidopsis_thaliana','brassica_rapa','glycine_max',
           'oryza_sativa','sorghum_bicolor','vitis_vinifera','zea_mays']

dict_l = []
r_path = '/home/yjko/data/plant/6_plant_numt_analysis/genome'

#numt list file name
nfn = 'numt_list_in_phylogenetic_tree_sort.txt'

def loading_db(sp_name):
    load_gff = open('%s/%s.gff3' %(r_path, sp_name),'r')
    load_aaseq = open('%s/%s_aa_h.seq' %(r_path, sp_name),'r')
    gff_r = load_gff.read()
    aaseq_r = load_aaseq.read()
    gffs = gff_r.strip().split('\n')
    header_l = []
    seq_l = []
    aaseqs = aaseq_r.split('\n>')
    for aaseq in aaseqs:
        info_l = aaseq.strip().replace('>','').split('\n')
        header = info_l[0]
        sequence = ''.join(info_l[1:])
        header_l.append(header)
        seq_l.append(sequence)
    aaseqs_dict = dict(zip(header_l,seq_l))
    return gffs, aaseqs_dict

handle = open('%s/%s' %(r_path, nfn),'r')
line = handle.readline()

added_sp_l = []

gff_db = []
aaseq_dic = {}

synteny_info = open('%s/numt_synteny_info.list' %(r_path),'w')
syn_seq = open('%s/numt_synteny_sequence.fa' %(r_path),'w')

while line != '':
    spline = line.split('\t')
    sp_n = spline[0]
    chr_n = spline[1]
    
    s_pos = int(spline[2])
    e_pos = int(spline[3])
    str_p = s_pos - fsl
    end_p = e_pos + fsl
    if s_pos < fsl:
        str_p = 1
    
    if sp_n not in added_sp_l:
        gff_db, aaseq_dic = loading_db(sp_n)
        added_sp_l.append(sp_n)
        
    bf_g = []
    bf_gs = []
    af_g = []
    af_gs = []
    for gff_s in gff_db:
        sgffs = gff_s.split('\t')
        chr_num = sgffs[0]
        categori = sgffs[2]
        g_stp = int(sgffs[3])
        g_edp = int(sgffs[4])
        
        if sp_n == 'sorghum_bicolor':
            chr_num = chr_num.replace('chromosome_','chr')
        elif sp_n != 'brassica_rapa':
            chr_num = 'chr%s' %(chr_num)
        
        if chr_num == chr_n:
            if str_p < g_edp:
                if g_stp < end_p:
                    #print str_p, end_p, g_stp, g_edp
                    g_id = ''
                    g_name = ''
                    desc = sgffs[8].split(';')
                    gname_list = []
                    
                    for x in desc:
                        if 'gene=' in x:
                            if 'EPl' not in x:
                                gn = x.split('=')[1].replace('transcript:','').replace('gene:','').upper()
                                if gn not in gname_list:
                                    gname_list.append(gn)
                        if 'transcript_id=' in x:
                            if 'EPl' not in x:
                                gn = x.split('=')[1].replace('transcript:','').replace('gene:','').upper()
                                if gn not in gname_list:
                                    gname_list.append(gn)

                        if 'Name=' in x:
                            if 'EPl' not in x:
                                gn = x.split('=')[1].replace('transcript:','').replace('gene:','').upper()
                                if gn not in gname_list:
                                    gname_list.append(gn)

                        if 'ID=' in x:
                            if 'EPl' not in x:
                                gn = x.split('=')[1].replace('transcript:','').replace('gene:','').upper()
                                if gn not in gname_list:
                                    gname_list.append(gn)
                    if len(gname_list) != 0:
                        if sp_n == 'oryza_sativa':
                            gname_list[0] = gname_list[0].replace('G','T')
                            if len(gname_list) ==2:
                                fst = gname_list[0]
                                snd = gname_list[1].replace('G','T')
                                if snd in fst:
                                    gname_list.remove(gname_list[1])
                        #print gname_list
                        #print '***************************************'
                        gene_id = ''
                        gene_name = ''
                        if len(gname_list) == 2:
                            gene_id = gname_list[0]
                            gene_name = gname_list[1]
                        else:
                            gene_id = gname_list[0]
                            gene_name = 'none'
                        g_sequence = aaseq_dic.get(gene_id)
                        
                        if g_edp <= s_pos:
                            bf_g.append(gene_id+'~'+gene_name)
                            bf_gs.append(g_sequence)
                        elif e_pos <= g_stp:
                            af_g.append(gene_id+'~'+gene_name)
                            af_gs.append(g_sequence)
                        else:
                            pass
    stl = line.strip().replace('\t','_')
    bfgl = '\t'.join(bf_g)
    afgl = '\t'.join(af_g)
    gene_synteny = '%s\t|----NUMT----|\t%s' %(bfgl, afgl)
    synteny_out_string = '%-45s\t%s\n' %(stl, gene_synteny)
    synteny_info.write(synteny_out_string)
    
    for num1 in range(len(bf_g)):
        head1 = bf_g[num1]
        seq1 = bf_gs[num1]
        out1_str = '>%s\n%s\n' %(head1, seq1)
        syn_seq.write(out1_str)
    syn_seq.write('\nNUMT ' + line + '\n')
    for num2 in range(len(af_g)):
        head2 = af_g[num2]
        seq2 = af_gs[num2]
        out2_str = '>%s\n%s\n' %(head2, seq2)
        syn_seq.write(out2_str)
    syn_seq.write('\n\n')
    line = handle.readline()
    
syn_seq.close()    
synteny_info.close()
