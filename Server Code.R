# Code for server

# protein_brain, protein_liver, secreted_protein_brain, secreted_protein_liver
# notsecreted_protein_brain, notsecreted_protein_liver

### Ssec calculations

target_liver = protein_liver

brain.liv = bicorAndPvalue(secreted_protein_brain, target_liver, use = 'pairwise.complete.obs')
# brain.liv$p

scores = rowSums(-log(brain.liv$p))

sec_Ssec_brain = data.frame(Gene_symbol = names(scores), score = scores) %>%
  #filter(Gene_symbol %in% sec_annot$ensembl_gene_id) %>%
  mutate(Ssec = score / length(colnames(target_liver))) %>%
  dplyr::select(-score) %>%
  arrange(desc(Ssec))


write.csv(sec_Ssec_brain, "Sec-Full-Ssec_brain.csv")
## LIVER ##

target_brain = protein_brain

liv.brain = bicorAndPvalue(secreted_protein_liver, target_brain, use = 'pairwise.complete.obs')
# brain.liv$p

scores = rowSums(-log(liv.brain$p))

sec_Ssec_liver = data.frame(Gene_symbol = names(scores), score = scores) %>%
  #filter(Gene_symbol %in% sec_annot$ensembl_gene_id) %>%
  mutate(Ssec = score / length(colnames(target_brain))) %>%
  dplyr::select(-score) %>%
  arrange(desc(Ssec))

write.csv(sec_Ssec_liver, "Sec-Full-Ssec_liver.csv")

#### Not Secreted ####

## BRAIN ##

brain.liv.c = bicorAndPvalue(notsecreted_protein_brain, target_liver, use = 'pairwise.complete.obs')
# brain.liv$p

scores = rowSums(-log(brain.liv.c$p))

nsec_Ssec_brain = data.frame(Gene_symbol = names(scores), score = scores) %>%
  #filter(Gene_symbol %in% sec_annot$ensembl_gene_id) %>%
  mutate(Ssec = score / length(colnames(target_liver))) %>%
  dplyr::select(-score) %>%
  arrange(desc(Ssec))

write.csv(nsec_Ssec_brain, "NSec-Full-Ssec_brain.csv")

## LIVER ##

liv.brain.c = bicorAndPvalue(notsecreted_protein_liver, target_brain, use = 'pairwise.complete.obs')
# brain.liv$p

scores = rowSums(-log(liv.brain.c$p))

nsec_Ssec_liver = data.frame(Gene_symbol = names(scores), score = scores) %>%
  #filter(Gene_symbol %in% sec_annot$ensembl_gene_id) %>%
  mutate(Ssec = score / length(colnames(target_brain))) %>%
  dplyr::select(-score) %>%
  arrange(desc(Ssec))

write.csv(nsec_Ssec_liver, "NSec-Full-Ssec_liver.csv")
