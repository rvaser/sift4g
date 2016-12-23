##Input options

The query file can have multiple protein sequences (in fasta format).

The substitution file(s) must have the same name as the ID in the protein file. Using the option --subst, pass in the **directory** containing the substitution file(s).  Based on the protein ID in the query file, SIFT 4G will look for the corresponding substitution file of the same name.

##Output

The predictions will be in &lt;Protein ID&gt;.SIFTprediction

To test SIFT 4G, go to the parent directory and type:

> ./bin/sift4g -q ./test_files/query.fasta --subst ./test_files/ -d ./test_files/sample_protein_database.fa

To get all 20 amino acid predictions for every position, don't pass in a substitution file:

> ./bin/sift4g -q ./test_files/query.fasta -d ./test_files/sample_protein_database.fa

Output (LACI_ECOLI.SIFTprediction and PURR_SALTY.SIFTprediction) will be in the directory where sift4g was run (unless specified otherwise by --out)                                                            
