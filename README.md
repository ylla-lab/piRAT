# piRNAnnotation
To run the code you have to (while being in the directory ../piRAT):<br>
```
pip3 install -r requirements.txt
```
And then<br>
```
pip3 install .
```
After this you can use this module <br>
```
pirat -p <directory_with_bam_files> -m primary
```
Args:<br>
```
  -h --help
  -p --path             Path to the directory with data, e.g data/ 
  -o --output-path      Path to the output directory, e.g output_data/
  -k --minreads         MinReads parameter 
  -e --eps              Eps parameter 
  -r --range_of_size    Range of size of piRNAs, e.g 26,32 
  -t --threads          Number of threads to use for clustering 
  -s --testing          USED FOR TESTING - Define number of scaffolds to perform clustering on them 
  -m --module           Which module to run: 'primary' for clustering algorithm, 'secondary' for ping-pong piRNA
                        detection, 'both' for both 
  -v --variation_threshold 
                        Value of nt in each direction of core read for cleansing loaded reads 
  -a                    If used, pirat will run automatically without asking user for validation of found parameters
  -d                    If used, graphical representation of found clusters will be generated
  --plot_iter           Number of plots generated per iteration. More == faster generation == more RAM use
  ```
Example use:<br>
```
pirat -p /data/ -k 8 -e 900 -r 28,29 -t 16 -d 
```