In this document I will explain how to use the software contained in this project. 

The easiest way to run the project is to execute the script named 'main.py'.  The defaults
that are provided are sufficient to run the project with the inputs provided.
Each accompanying script can be run on its own as well.  main.py does not allow for the intermediate
files or folders to be renamed by the user.  Currently the user cannot select the point that analysis
begins at through main.py, but this functionality will be added later.

./main.py -h
usage: main.py [-h] [-i FILE] [-I FOLDER] [-o FOLDER] [-f FILE] [-n INT]
               [-m INT] [-g INT] [-e FLOAT]

The purpose of this script is to run the full software suite that we have
developed to study operons using as few inputs as possible. This will
facilitate the ease of use as much as possible.

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --infile FILE
                        Input file for the operon query step of the pipeline.
  -I FOLDER, --infolder FOLDER
                        Folder containing all genbank files for use by the
                        program.
  -o FOLDER, --outfolder FOLDER
                        Folder where results will be stored.
  -f FILE, --filter FILE
                        File restrictiong which accession numbers this script
                        will process. If no file is provided, filtering is not
                        performed.
  -n INT, --num_proc INT
                        Number of processors that you want this script to run
                        on. The default is every CPU that the system has.
  -m INT, --min_genes INT
                        Minum number of genes that an operon must contain
                        before it can be considered for further analysis. The
                        default is 5 because that is what we are currently
                        using in the study.
  -g INT, --max_gap INT
                        Size in nucleotides of the maximum gap allowed between
                        genes to be considered neighboring. The default is
                        500.
  -e FLOAT, --eval FLOAT
                        eval for the BLAST search.
