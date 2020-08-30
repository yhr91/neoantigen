library("biomaRt")
library("stringr")
library("seqinr")
library("Biostrings")
ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl",host = "uswest.ensembl.org")
#ensembl = useEnsembl("ensembl","hsapiens_gene_ensembl", mirror = "uswest")

library(rtracklayer)
setwd("C:\\Users\\trabb\\Documents")
#gtf file was filtered using bedtools to select only the exons coding for protein coding genes
gtf = rtracklayer::import('GRCh38.99_codingExons.gtf')
gtf_df=as.data.frame(gtf)
gtf_df = gtf_df[which(gtf_df$transcript_biotype == "protein_coding"),]

bed = read.table('Genome_Canonical.bed', header = FALSE, sep = '\t')
bed_df = as.data.frame(bed)
bed_df$V5 = as.character(bed_df$V5)

#removes the transcript version from the bed file for easier matching with the GTF
#need to do this so as to find matches
for (i in 1:length(bed_df$V5)){
  bed_df$V5[i] = substr(bed_df$V5[i],1,15)
}

#find all the canonical transcripts
matches = gtf_df$transcript_id %in% bed_df$V5
matchGTF = gtf_df[matches,]
unq_geneIDs = unique(matchGTF$gene_id)

#these columns are type int and need to be coerced to character so that they can be edited
matchGTF$type = as.character(matchGTF$type)
matchGTF$seqnames = as.character(matchGTF$seqnames)
matchGTF$transcript_id = as.character(matchGTF$transcript_id)
matchGTF$transcript_version = as.character(matchGTF$transcript_version)

#since the 'type' column is 'exon' for all of them, change it to a useful fasta header for use in bedtools and later on in this script
for (i in 1:length(matchGTF$seqnames)){
  gene_id = paste(c("gene_id:",matchGTF$gene_id[i]),collapse = "")
  trans_id = paste(c("transcript_id:",matchGTF$transcript_id[i]),collapse = "")
  #exon_no = paste(c("exon_number:",matchGTF$exon_number[i]),collapse = "")
  gene_name = paste(c("gene_name:",matchGTF$gene_name[i]),collapse = "")
  #name = paste(c(gene_id,trans_id,exon_no,gene_name),collapse = "|")
  name = paste(c(gene_id,trans_id,gene_name),collapse = "|")
  matchGTF$type[i] = name
}


##################################################
#Main loop
#the main loop of this script - gets all possible frameshifts for each PAM for each gene
for(i in 1:length(unq_geneIDs)){
  #get all exons for a specific gene
  geneI = matchGTF[matchGTF$gene_id %in% unq_geneIDs[i],]
  
  #retrieve the exon seqs for a particular gene. need to do this in a for loop b/c providing a list of IDs returns the exon sequences in a random order
  exon_seqs = c()
  for(j in 1:length(geneI$exon_id)){
    exon_id = geneI$exon_id[j]
    exon_seq = biomaRt::getSequence(type = "ensembl_exon_id",id = exon_id,seqType = "gene_exon",mart = ensembl)
    exon_seqs = c(exon_seqs,exon_seq$gene_exon)
  }
  
  noTerminalExon = DNAString(paste(exon_seqs[1:(length(exon_seqs)-1)],collapse = ""))
  #get the cDNA sequence for the same gene
  trans_id = paste(c(geneI$transcript_id[1],geneI$transcript_version[1]),collapse = ".")
  cds = biomaRt::getSequence(id = trans_id, type = "ensembl_transcript_id_version",seqType = "coding",mart = ensembl)
  cdna = biomaRt::getSequence(id = trans_id, type = "ensembl_transcript_id_version",seqType = "cdna",mart = ensembl)
  #gene_flank = biomaRt::getSequence(id = trans_id,type = "ensembl_transcript_id_version",seqType = "transcript_flank", upstream = 20,mart = ensembl)
  
  cds_dna = DNAString(x = cds$coding)
  cdna_dna = DNAString(x = cdna$cdna)
  #gene_flank_dna = DNAString(x = gene_flank$transcript_flank)
  #cdna_dna = xscat(gene_flank_dna,cdna_dna) #adding the upstream sequence for determing sgRNA
  
  #when the sgRNA is at the 5' or 3' ends of the coding region, we may need the UTRs from the cDNA to determing the sgRNA sequence
  offset = start(matchPattern(pattern = cds_dna[1:21],subject = cdna_dna))
  
  forwardStrandPAMs = matchPattern(pattern = "NGG",subject = cds_dna,fixed = FALSE)
  #if the PAM is any closer to the 5' end, then an indel will disrupt the start codon. This could eliminate any translation product or lead
  #to the use of a downstream start codon. Either way, considering this biology complicates the analysis. For simplicity, we're
  #only considering the guides that can lead to neopeptides using the canonical start codon.
  forwardStrandPAMs = forwardStrandPAMs[which(start(forwardStrandPAMs)>8)]
  
  reverseStrandPAMs = matchPattern(pattern = "CCN",subject = cds_dna,fixed = FALSE)
  #if the PAM on the reverse strand is the codon next to the stop codon, then the edit will only impact the 3'UTR and not the coding sequence,
  #meaning that future code will not work. Filter those that are more than a codon away. Additionally, there are some genes that don't have
  #(annotated) 3'UTRs (i.e. cDNA ends at stop), so code will fail unless you restrict PAMs to >= 8bp away from end
  reverseStrandPAMs = reverseStrandPAMs[which(start(reverseStrandPAMs) <= (length(cds_dna)-8))]
  
  #now that we have all the PAMs in the coding sequence, we can concatenate the 3'UTR on the end of the coding sequence
  #this is needed if any of the guides lead to a frameshift into the UTR
  #some genes don't have 3'UTR (or 5'UTR for that matter), so need to take that into consideration
  if(length(cdna_dna[offset:length(cdna_dna)]) == length(cds_dna)){
    #there's no 3'UTR, so don't update the coding sequence DNA variable
  } else{
    cds_dna = xscat(cds_dna,cdna_dna[(length(cds_dna)+offset):length(cdna_dna)])
  }
  
  
  #first all the PAMs on the forward strand
  #####
  #for all -1/-2 indels, the de novo reading frame begins at the codon before the PAM. This acts as an offset to get only the de novo
  #protein sequence for neoantigen prediction
  shift = c(5,3,4)
  shift_rev = c(4,6,5)
  #go through each of the PAMs
  if(!isEmpty(forwardStrandPAMs)){
    for(k in 1:length(forwardStrandPAMs)){
      temp = forwardStrandPAMs[k]
      minusOne = xscat(cds_dna[1:(start(temp)-5)],cds_dna[(start(temp)-3):length(cds_dna)])
      minusTwo = xscat(cds_dna[1:(start(temp)-6)],cds_dna[(start(temp)-3):length(cds_dna)])
      
      #gets the first stop codon in the -1 frameshift
      firstStop = min(c(which(codons(minusOne) == DNAString("TGA")), which(codons(minusOne) == DNAString("TAA")), which(codons(minusOne) == DNAString("TAG"))))
      
      if(is.infinite(firstStop)){
        #there's no stop codon in the transcript -> non-stop decay! shouldn't see a protein in this scenario
      } else{
        start = max(start(temp)-21+offset, 1)
        end = start(temp)-2+offset
        m1start = start(temp)-1
        m1pep = Biostrings::translate(minusOne[(m1start - shift[(m1start %% 3)+1]):(firstStop*3-3)])
        
        if(length(m1pep)<8){
          #if the window isn't present in the last two exons, then it won't lead to neoantigen, test minusTwo
        } else if(length(exon_seqs)==1){
          sgRNA = paste(c("sgRNA",as.character(cdna_dna[start:end])),collapse = ":")
          name1 = paste(c(geneI$type[1],"m1fs",sgRNA,"strand:+"),collapse = "|")
          write.fasta(sequences = list(m1pep),names = c(name1),open = "a",file.out = "frameshiftPeptides.fasta")
        } else if((length(noTerminalExon) - (firstStop*3-2+offset)) < 50){
          #it's in the last two and it's within 50bp of the exon-exon junction, so should bypass NMD
          sgRNA = paste(c("sgRNA",as.character(cdna_dna[start:end])),collapse = ":")
          name1 = paste(c(geneI$type[1],"m1fs",sgRNA,"strand:+"),collapse = "|")
          write.fasta(sequences = list(m1pep),names = c(name1),open = "a",file.out = "frameshiftPeptides.fasta")
        }
      }
      
      
      
      
      #gets the first stop codon in the -2 frameshift
      secondStop = min(c(which(codons(minusTwo) == DNAString("TGA")), which(codons(minusTwo) == DNAString("TAA")), which(codons(minusTwo) == DNAString("TAG"))))
      if(is.infinite(secondStop)){
        #non-stop decay, skip
      } else{
        start = max(start(temp)-21+offset,1)
        end = start(temp)-2+offset
        m2start = start(temp)-2
        m2pep = Biostrings::translate(minusTwo[(m2start - shift[(m2start %% 3)+1]):(secondStop*3-3)])
        
        if(length(m2pep)<8){
          #if the window isn't present in the last two exons, then it won't lead to neoantigen, test minusTwo
        } else if(length(exon_seqs)==1){
          #if the gene is a single-exon gene
          sgRNA = paste(c("sgRNA",as.character(cdna_dna[start:end])),collapse = ":")
          name2 = paste(c(geneI$type[1],"m2fs",sgRNA,"strand:+"),collapse = "|")
          write.fasta(sequences = list(m2pep),names = c(name2),open = "a",file.out = "frameshiftPeptides.fasta")
        } else if((length(noTerminalExon) - (secondStop*3-2+offset)) < 50){
          #it's in the last two and it's within 50bp of the exon-exon junction, so should bypass NMD
          sgRNA = paste(c("sgRNA",as.character(cdna_dna[start:end])),collapse = ":")
          name2 = paste(c(geneI$type[1],"m2fs",sgRNA,"strand:+"),collapse = "|")
          write.fasta(sequences = list(m2pep),names = c(name2),open = "a",file.out = "frameshiftPeptides.fasta")
        }
        
        
      }
    }
  }
  
  
  
  #now need to do all of this again for the PAMs on the reverse strand
  #####
  if(!isEmpty(reverseStrandPAMs)){
    for(l in 1:length(reverseStrandPAMs)){
      temp = reverseStrandPAMs[l]
      #could make indels on the reverse complement, then take the reverse complement to get them back into the normal reading frame,
      #but easier to search for CCN and work with that
      minusOne = xscat(cds_dna[1:(start(temp)+5)],cds_dna[(start(temp)+7):length(cds_dna)])
      minusTwo = xscat(cds_dna[1:(start(temp)+5)],cds_dna[(start(temp)+8):length(cds_dna)])
      
      firstStop = min(c(which(codons(minusOne) == DNAString("TGA")), which(codons(minusOne) == DNAString("TAA")), which(codons(minusOne) == DNAString("TAG"))))
      if(is.infinite(firstStop)){
        #non-stop decay, skip
      } else{
        start = start(temp)+2+offset
        end = min(start(temp)+21+offset,length(minusOne)+offset)
        m1start = start(temp)
        m1pep = Biostrings::translate(minusOne[(m1start + shift_rev[(m1start %% 3)+1]):(firstStop*3-3)])
        
        if(length(m1pep)<8){
          #will only search for neoantigens from length 8-12aa, so don't include anything smaller than 8
        } else if(length(exon_seqs)==1){
          #if the gene is a single-exon gene
          sgRNA = paste(c("sgRNA",as.character(reverseComplement(cdna_dna[start:end]))),collapse = ":")
          name1 = paste(c(geneI$type[1],"m1fs",sgRNA,"strand:-"),collapse = "|")
          write.fasta(sequences = list(m1pep),names = c(name1),open = "a",file.out = "frameshiftPeptides.fasta")
        } else if((length(noTerminalExon) - (firstStop*3-2+offset)) < 50){
          #need to take the reverse complement of the coding sequence DNA with adjusted offset to get sgRNA sequence
          sgRNA = paste(c("sgRNA",as.character(reverseComplement(cdna_dna[start:end]))),collapse = ":")
          name1 = paste(c(geneI$type[1],"m1fs",sgRNA,"strand:-"),collapse = "|")
          write.fasta(sequences = list(m1pep),names = c(name1),open = "a",file.out = "frameshiftPeptides.fasta")
        }
      }
      
      
      
      secondStop = min(c(which(codons(minusTwo) == DNAString("TGA")), which(codons(minusTwo) == DNAString("TAA")), which(codons(minusTwo) == DNAString("TAG"))))
      if(is.infinite(secondStop)){
        #result is non-stop decay, skip
      } else{
        start = start(temp)+2+offset
        end = min(start(temp)+21+offset,length(minusTwo)+offset)
        m2start = start(temp)
        m2pep = Biostrings::translate(minusTwo[(m2start + shift_rev[(m2start %% 3)+1]):(secondStop*3-3)])
        
        if(length(m2pep)<8){
          #
        } else if(length(exon_seqs)==1){
          sgRNA = paste(c("sgRNA",as.character(reverseComplement(cdna_dna[start:end]))),collapse = ":")
          name2 = paste(c(geneI$type[1],"m2fs",sgRNA,"strand:-"),collapse = "|")
          write.fasta(sequences = list(m2pep),names = c(name2),open = "a",file.out = "frameshiftPeptides.fasta")
        } else if((length(noTerminalExon) - (secondStop*3-2+offset)) < 50){
          sgRNA = paste(c("sgRNA",as.character(reverseComplement(cdna_dna[start:end]))),collapse = ":")
          name2 = paste(c(geneI$type[1],"m2fs",sgRNA,"strand:-"),collapse = "|")
          write.fasta(sequences = list(m2pep),names = c(name2),open = "a",file.out = "frameshiftPeptides.fasta")
        }
      }
      
      
      
    }
  }
  
}

#BED file changes
#ENST00000379099.3 -> ENST00000379116.10
#remove ENSG00000284956 (incomplete CDS)
#remove ENSG00000236737 (no ATG start)
#remove ENSG00000257390 (no ATG start)
#remove ENSG00000275674 (no ATG start)
#remove ENSG00000261832 (incomplete CDS)
#skip ENSG00000284741 (code fails on Sherlock, use local results)
#reverse strand sgRNA seqs might be wrong for some genes in fasta 1,2-4,6-10,12,13; check later
#remove ENSG00000285625 (incomplete CDS)
#remove ENSG00000288000 (incomplete CDS)
#skip ENSG00000285133 for now
#remove ENSG00000257921
#remove ENSG00000286143 (incomplete CDS)
#remove ENSG00000271810 (incomplete CDS)
#remove ENSG00000273154 (incomplete CDS)
#remove ENSG00000274810 (incomplete CDS)
#remove ENSG00000257767 (incomplete CDS)
#remove ENSG00000278646 (incomplete CDS)
#remove ENSG00000249141 (incomplete CDS)
#remove ENSG00000283149 (incomplete CDS)
#remove ENSG00000250424 (incomplete CDS)
#remove ENSG00000259132 (incomplete CDS)
#remove ENSG00000285947 (incomplete CDS)
#remove ENSG00000259371 (incomplete CDS)
#remove ENSG00000161939 (incomplete CDS)
#remove ENSG00000273047 (incomplete CDS)
#replace ENST00000586691.1 with ENST00000317702.10
#remove ENSG00000270249 (incomplete CDS)
#remove ENSG00000205592 (incomplete CDS)
#remove ENSG00000272297 (incomplete CDS)
#replace ENST00000639600.1 with ENST00000640252.2
#remove ENSG00000261740 (incomplete CDS)
#remove ENSG00000258989 (incomplete CDS)
#skip ENSG00000114270 (BioMart doesn't work for this gene for whatever reason)
#replace ENST00000273980.9 with ENST00000394708.7
#skip ENSG00000054654 (can't retrieve exon 13)
#remove ENSG00000283528 (incomplete CDS)
#remove ENSG00000267740 (incomplete CDS)
#remove ENSG00000268861 (incomplete CDS)
#ENSG00000267952
#ENSG00000269711
#ENSG00000268870
#ENSG00000269590
#ENSG00000269095
#ENSG00000269035
#ENSG00000283201
#ENSG00000266953
#ENSG00000267360
#ENSG00000267748
#ENSG00000268083
#ENSG00000269547
#ENSG00000267881
#ENSG00000268643
#ENSG00000268361
#ENSG00000267173
#ENSG00000268434
#ENSG00000268465
#ENSG00000267335
#ENSG00000268133
#ENSG00000269533
#ENSG00000268533
#ENSG00000269026
#ENSG00000276612
#ENSG00000280433
#ENSG00000173366
#replace ENST00000431216.5 with ENST00000673807.1
#ENSG00000249624
#replace ENST00000479817 with ENST00000252971.11
#replace ENST00000429029 with ENST00000262177.9	
#change ENST00000389394.7 to ENST00000389394.8
#ENSG00000261793
#replace ENST00000342992 with ENST00000589042.5
#ENSG00000163602
#ENSG00000180264
#need to run ENSG00000155657 (change ENST00000342992.10 to ENST00000589042.5), ENSG00000183091 later
#ENSG00000142539
#replace ENST00000409184.7 with ENST00000375458.6
#replace ENST00000541798.1 with ENST00000281938.6
#do ENSG00000060718, ENSG00000143341, (1.bed)
#do ENSG00000065457 later (22)
#do ENSG00000072657 later (18)
#replace ENST00000613249.4 with ENST00000356725.9
#ENSG00000282246
#ENSG00000283496
#ENSG00000111780
#skip ENSG00000131018 (9) ENSG00000154358 (1)
#ENSG00000131152
#skip ENSG00000162825 (1)
#replace ENST00000451047.1 with ENST00000274979.12
#replace ENST00000595048.5 with ENST00000435989.7
#replace ENST00000596929 with ENST00000269829.5
#replace ENST00000367614.5 with ENST00000444136.6