#%%

def classify_protein_splice_fsm(row):
    if row.pr_cat=='full-splice_match' and row.pr_subcat=='multi-exon':
        if row.pr_nterm_diff==0 and row.pr_cterm_diff==0:
            return 'pFSM,,,'
        elif (row.pr_nterm_diff!=0 or row.pr_cterm_diff!=0) and row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNIC,combo_nterm_cterm,,is alt nterm or cterm?'
        elif row.tx_cat=='incomplete-splice_match' and row.pr_ntermhang<0 and row.pr_ctermhang==0:
            return 'pISM,cand_ntrunc,,assume nterm novel'
        elif row.tx_cat=='full-splice_match' and row.tx_5hang<0 and row.pr_ntermhang>0:
            return 'pISM,cand_ntrunc,,assume nterm novel'
        elif (row.tx_cat=='novel_in_catalog' or row.tx_cat=='novel_not_in_catalog') and row.pr_ntermhang<0 and row.pr_ctermhang==0:
            return 'pNIC,alt_nterm,,assume nterm novel'
        elif (row.tx_cat=='full-splice_match' or row.tx_cat=='incomplete-splice_match') and row.tx_5hang>0 and row.pr_ntermhang<0 and row.pr_ctermhang==0:
            return 'pNIC,alt_nterm,,assume nterm novel'
        elif row.pr_nterm_diff==0 and row.pr_cterm_diff!=0:
            return 'pISM,cand_ctrunc,,'
        elif row.tx_cat=='full-splice_match' and row.pr_cterm_diff==0 and row.tx_5hang<0 and row.pr_ntermhang<0:
            return 'pISM,cand_ntrunc,,'
        elif row.pr_ntermhang<0 and row.pr_ctermhang<0:
            return 'unicorn?,,,'
        elif row.pr_ntermhang>0 and row.pr_ctermhang==0:
            return 'pNIC,alt_nterm,,'
        elif row.pr_nterm_gene_diff!=0 and row.pr_ntermhang>0:
            return 'pNNC,novel_nterm_known_splice_known_cterm,,'
        elif row.pr_nterm_gene_diff==0 and row.pr_ctermhang>0:
            return 'unicorn?,,,'
        elif row.pr_ntermhang>0 and row.pr_ctermhang>0:
            return 'unicorn?,,,'
        else:
            return 'orphan'
    return ''


def classify_protein_splice_ism(row):
    if row.pr_cat=='incomplete-splice_match' and 'mono-exon' not in row.pr_subcat:
        if row.pr_nterm_diff==0 and row.pr_cterm_diff!=0 and row.pr_ctermhang<0:
            return 'pISM,cand_ctrunc,,'
        elif row.pr_nterm_diff==0 and row.pr_cterm_diff!=0 and row.pr_cterm_gene_diff!=0 and row.pr_ctermhang>0:
            return 'pNNC,known_nterm_known_splice_novel_cterm,is_nmd?,'
        elif row.pr_nterm_diff==0 and row.pr_cterm_diff!=0 and row.pr_cterm_gene_diff==0 and row.pr_ctermhang>0:
            return 'pNIC,known_nterm_known_splice_alt_cterm,is_nmd?,'
        elif row.tx_cat=='incomplete-splice_match' and row.pr_nterm_diff!=0 and row.tx_5hang<0 and row.pr_ntermhang<0:
            return 'pISM,cand_ntrunc,,'
        elif row.tx_cat=='incomplete-splice_match' and row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.tx_5hang<=10 and row.pr_ntermhang<0:
            return 'pISM,cand_ntrunc,,'
        elif row.tx_cat=='full-splice_match' and row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_ntermhang<0:
            return 'pISM,cand_ntrunc,,'
        elif row.tx_cat=='novel_in_catalog' and row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_nterm_gene_diff!=0 and row.pr_ntermhang<0:
            return 'pNNC,novel_nterm_known_splice_known_cterm,,'
        elif row.tx_cat=='novel_not_in_catalog' and row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_nterm_gene_diff==0 and row.pr_ntermhang<0:
            return 'pNIC,alt_nterm_known_splice_known_cterm,,'
        elif row.tx_cat=='novel_not_in_catalog' and row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_nterm_gene_diff!=0 and row.pr_ntermhang<0:
            return 'pNNC,novel_nterm_known_splice_known_cterm,,'
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_nterm_gene_diff==0 and row.pr_ntermhang>0:
            return 'pNIC,alt_nterm_known_splice_known_cterm,,'
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff==0 and row.pr_nterm_gene_diff!=0 and row.pr_ntermhang>0:
            return 'pNNC,novel_nterm_known_splice_known_cterm,,'
        elif row.pr_nterm_diff!=0 and row.pr_cterm_diff!=0:
            return 'unicorn?,,,'
        else:
            return 'orphan_ism'
    return ''
    

def classify_protein_splice_nic(row):
    if row.pr_cat=='novel_in_catalog' and 'mono-exon' not in row.pr_subcat:
        if row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNIC,known_nterm_combo_splice_known_cterm,,'
        elif row.pr_nterm_gene_diff!=0 and row.pr_cterm_gene_diff==0:
            return 'pNNC,novel_nterm_combo_splice_known_cterm,,'
        elif row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,known_nterm_combo_splice_novel_cterm,is_nmd?,'
        elif row.pr_nterm_gene_diff!=0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,novel_nterm_combo_splice_novel_cterm,is_nmd?,'
        else:
            return 'orphan_nic'
    return ''


def classify_protein_splice_nnc(row):
    if row.pr_cat=='novel_not_in_catalog':
        if row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff==0:
            return 'pNNC,known_nterm_novel_splice_known_cterm,,'
        elif row.pr_nterm_gene_diff!=0 and row.pr_cterm_gene_diff==0:
            return "pNNC,novel_nterm_novel_splice_known_cterm,,could be FL or 5' deg product"
        elif row.pr_nterm_gene_diff==0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,known_nterm_novel_splice_novel_cterm,is_nmd?,could be FL or NMD product'
        elif row.pr_nterm_gene_diff!=0 and row.pr_cterm_gene_diff!=0:
            return 'pNNC,novel_nterm_novel_splice_novel_cterm,is_nmd?,,'
        else:
            return 'orphan_nnc'
    return ''


=IF(OR(E2="mono-exon",E2="mono-exon_by_intron_retention"),
IF(AND(C2="full-splice_match",J2=0,K2=0),"pFSM,mono-exon,,",
IF(C3="intergenic","intergenic,mono-exon,,",
IF(L2=0, "trunc/altnterm","orphan_monoexon"))),"")
def classify_protein_splice_monoexon(row):
    if row.pr_subcat=='mono-exon' or row.pr_subcat=='mono-exon_by_intron_retention':
        if row.pr_cat=='full-splice_match' and row.pr_nterm_diff==0 and row.pr_cterm_diff==0:
            return 'pFSM,mono-exon,,'
        elif row.pr_cat=='intergenic':
            return 'intergenic,mono-exon,,'
        elif row.pr_nterm_gene_diff==0:
            return 'trunc/altnterm'
        else:
            return 'orphan_monoexon'
    return ''

def classify_protein_splice_misc(row):
    if row.pr_subcat=='multi-exon':
        if row.pr_cat=='intergenic':
            return 'intergenic,multi-exon,,'
        elif row.pr_cat=='fusion':
            return 'fusion,multi-exon,,'
    return ''


def classify_protein(row):
    fsm_classiciation = classify_protein_splice_fsm(row)
    ism_classificaiton = classify_protein_splice_ism(row)
    nic_classification = classify_protein_splice_nic(row)
    nnc_classification = classify_protein_splice_nnc(row)
    mono_classification = classify_protein_splice_monoexon(row)
    misc_classification = classify_protein_splice_misc(row)
    return (
        fsm_classiciation+
        ism_classificaiton+
        nic_classification+
        nnc_classification+
        mono_classification+
        misc_classification
    )
    

# %%
