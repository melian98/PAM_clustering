#Make a list of the required libraries
libraries <- c("rentrez", "BiocManager", "dplyr", "tidyverse", "factoextra", "hopkins", "fpc", "cluster")

#initialize needed variables
i = 1
less_2000_TNF = 0
less_3000_TNF = 0
less_4000_TNF = 0
less_5000_TNF = 0
less_2000_P53 = 0
less_3000_P53 = 0
less_4000_P53 = 0
less_5000_P53 = 0

#install packages if they're not currently installed
install.packages(setdiff(libraries, rownames(installed.packages())))
BiocManager::install("Biostrings")

#load the needed librariers
lapply(libraries, library, character.only = TRUE)
library(Biostrings)

#set a goal for the silhouette value (between 0 and 1, higher is better). A score of 0.5 or greater is considered good and increasing this value may greatly increase computation time
goal = 0.6

#initial silhouette value
silhouette_value = 0

#if the user does not have TNF fetch and P53 fetch fasta files already prepared, they will be downloaded
if (!("TNF_fetch.fasta" %in% list.files() & "P53_fetch.fasta" %in% list.files())) {
  
  #our search of the NCBI database will query the nucleotide database for TNF genes between 1000 and 10 000 base pairs   long
  NCBI_TNF <- entrez_search(db = "nuccore", term = "TNF[Gene Name] OR tnfa[Gene Name] OR tnfb[Gene Name] AND 1000:10000[SLEN]")
  NCBI_P53 <- entrez_search(db = "nuccore", term = "(P53[Gene Name]) OR TP53[Gene Name] AND 1000:5000[SLEN] ") 
  
  #save the retmax count in order to repeat the search with this new count
  new_retmax_TNF <- NCBI_TNF$count
  new_retmax_P53 <- NCBI_P53$count
  
  #perform the search again with the new retmax value
  NCBI_TNF <- entrez_search(db = "nuccore", term = "TNF[Gene Name] OR tnfa[Gene Name] OR tnfb[Gene Name] AND 1000:10000[SLEN]", retmax = new_retmax_TNF, use_history = TRUE)
  
  NCBI_P53 <- entrez_search(db = "nuccore", term = "(P53[Gene Name]) OR TP53[Gene Name] AND 1000:5000[SLEN]", retmax = new_retmax_P53, use_history = TRUE)
  
  #retrieve the results from the web search in order to save them to a fasta file
  TNF_fetch <- entrez_fetch(db = "nuccore", rettype = "fasta", web_history = NCBI_TNF$web_history)
  P53_fetch <- entrez_fetch(db = "nuccore", rettype = "fasta", web_history = NCBI_P53$web_history)
  
  #write the fetch variables to fasta files
  write(TNF_fetch, "TNF_fetch.fasta", sep = "\n")
  write(P53_fetch, "P53_fetch.fasta", sep = "\n")
  
}

#retrieve the fasta files
TNF_stringset <- readDNAStringSet("TNF_fetch.fasta")
P53_stringset <- readDNAStringSet("P53_fetch.fasta")

#paste the sequences from the fasta file onto a dataframe
dataframe_TNF <- data.frame(sequence = paste(TNF_stringset))
dataframe_P53 <- data.frame(sequence = paste(P53_stringset))

#check the length of the nucleotide sequences in TNF to look for any outliers
for (k in 1:length(dataframe_TNF$sequence)){

  #let seq_length represent the number of nucleotides in each sequence
  seq_length <- nchar(dataframe_TNF$sequence[k])
  
  #increment the appropriate counter depending on the length of the nucleotide sequence
  if (seq_length < 2000) {less_2000_TNF = less_2000_TNF + 1}
  else if (seq_length < 3000) {less_3000_TNF = less_3000_TNF + 1}
  else if (seq_length < 4000) {less_4000_TNF = less_4000_TNF + 1}
  else {less_5000_TNF = less_5000_TNF + 1}
}

#save the TNF counters into a list
length_TNF <- c(less_2000_TNF, less_3000_TNF, less_4000_TNF, less_5000_TNF)

#check the length of the nucleotide sequences in P53 to look for any outliers
for (k in 1:length(dataframe_P53$sequence)){
  
  #let seq_length represent the number of nucleotides in each sequence
  seq_length <- nchar(dataframe_P53$sequence[k])
  
  #increment the appropriate counter depending on the length of the nucleotide sequence
  if (seq_length < 2000) {less_2000_P53 = less_2000_P53 + 1}
  else if (seq_length < 3000) {less_3000_P53 = less_3000_P53 + 1}
  else if (seq_length < 4000) {less_4000_P53 = less_4000_P53 + 1}
  else {less_5000_P53 = less_5000_P53 + 1}
}

#save the P53 counters into a list
length_P53 <- c(less_2000_P53, less_3000_P53, less_4000_P53, less_5000_P53)

#save the two lists into a dataframe for plotting
counts <- data.frame(length_TNF, length_P53)

#transpose the counts dataframe
counts <- t(counts)

#convert the counts dataframe into a matrix so it can be used in the barplot
counts <- data.matrix(counts)

#create a barplot to visualize any outliers in nucleotide length
barplot(counts, main = "frequency of different sequence length", xlab = "sequence length", names.arg = c("less_2000", "less_3000", "less_4000", "less_5000"), ylab = "sequence length (nt)", col = c("darkblue", "red"), beside = TRUE)

#create a legend in the barplot for clarity
legend("topright", c("length_P53", "length_TNF"), fill = c("darkblue", "red"))

#remove all unnecessary variables used for NCBI search and barplot creation
rm(less_2000_TNF, less_3000_TNF, less_4000_TNF, less_5000_TNF, less_2000_P53, less_3000_P53, less_4000_P53, less_5000_P53, NCBI_P53, new_retmax_P53, P53_fetch, P53_stringset, NCBI_TNF, new_retmax_TNF, TNF_fetch, TNF_stringset, libraries, counts, k, length_TNF, length_P53, seq_length)

#set the seed before continuing to ensure reproducible results
set.seed(100)

#take a sample of the larger dataframe equal in size to the smaller dataframe in order to have an equal number of results for both dataframes
dataframe_P53 <- sample_n(dataframe_P53, size = length(dataframe_TNF$sequence))

#combine both dataframes into a single large dataframe
combined_dataframe_original <- rbind(dataframe_TNF, dataframe_P53)

#remove the previous two dataframes as they are no longer needed
rm (dataframe_TNF, dataframe_P53)

#remove any whitespaces that may exist in the new dataframe
combined_dataframe_original$sequence <- gsub(" ", "", combined_dataframe_original$sequence)

#scramble the combined dataframe to ensure that the results from both genes are no longer dinstinguishable from each other
combined_dataframe_original <- sample_frac(combined_dataframe_original, 1L)

#convert the sequence information to an object of type DNAStringset
combined_dataframe_original$sequence <- DNAStringSet(combined_dataframe_original$sequence)

#This loop will ensure that the silhouette value reaches the users goal
while (silhouette_value < goal) {
  #combined_dataframe takes the value from the original dataframe before the loop began
  combined_dataframe <- combined_dataframe_original
  #combined dataframe will join columns to the dataframe with an oligonucleotide of length i
  combined_dataframe <- cbind(combined_dataframe, as.data.frame(oligonucleotideFrequency(combined_dataframe$sequence, width = i)))
  
  #Use the PAM clustering algorithm to cluster the columns with the count of the oligonucleotides
  clustering_pam <- pam(combined_dataframe[ , 2:ncol(combined_dataframe)], 2, metric = "euclidean", stand = FALSE)

  #retrieve the silhouette value from the clustering algorithm to check if it has reached the users goal
  silhouette_info <- clustering_pam$silinfo
  silhouette_value <- silhouette_info$avg.width
  
  #increment i by 1 and repeat the loop if the goal has not been met
  i = i+1
}

#create a cluster plot to show the cluster groups made from the PAM algorithm
fviz_cluster(clustering_pam, data = combined_dataframe[ , 2:ncol(combined_dataframe)], main = "Gene Clusters") 

#create a dissimilarity matrix to visualize the disorder from the PAM clustering compared to the disorder of a random sample
fviz_dist(dist(combined_dataframe[ , 2:ncol(combined_dataframe)]), show_labels = FALSE)+ labs(title = "cluster dissimilarity")

#generate a random sample to compare the disorder
random <- replicate(6, sample (1:60, 40, rep=TRUE))

#plot another dissimilarity matrix for the random sample
fviz_dist(dist(random), show_labels = FALSE)+ labs(title = "random_data")

#it is clear that there is a more ordered pattern in the dissimilarity matrix produced from the clustering data indicating that there is order in the clustering

# a silhoutte plot is also generated to visually observe the silhouette index
fviz_silhouette(clustering_pam, palette = "jco", ggtheme = theme_classic())

#check the Calinski-Harabasz index to confirm the quality of the clustering.
calinhara <- calinhara(combined_dataframe[ , 2:ncol(combined_dataframe)], clustering = clustering_pam$cluster, cn = 2)
