# herd-immunity
Code for paper "Investigating vaccine-induced immunity and its effect in mitigating SARS-CoV-2 epidemics in China"

This is under review and available on MedrXiv: https://www.medrxiv.org/content/10.1101/2021.07.23.21261013v2

The herd immunity.rar includes all the materials to generate the most results in the above paper, while the rest simply results can be generated by adopting the provided codes.  

The files are described as follows:

1. The code directory includes two .cpp files (all-or-nothing model.cpp, and leaky.cpp), corresponding the two vaccine model used in the above paper, and one .sh file (sample.sh) which is an example of input parameters.

2. The input director includes all the data used in the paper, namely  
china_contact-matrix: the 16\times 16 contact matrix of China;  
contraindication: the proportion of contraindication and pregnant women by age groups;  
dose: the historical vaccine administration in China;  
population: the shanghai population by age groups;  
susceptibility: the age-specific susceptibility;  
shanghai-contact-matrix: includes 200 bootstrapped 16\times 16 contact matrices of Shanghai.  

If any bugs or issues with the code are identified please do get in touch!




