# Relion_star_parser
A simple Relion .star file parser which allows Relion star file plotting, modification and conversion.

It read star file and converts this to the Pandas DataFrame what simplifies operations on the data, and finally saves star file.

# Requirements
To run the script some external libraries are required:

Python 3+ with:
+ matplotlib (plotting selected rows)
+ gemmi (provides star and cif file parser)
+ numpy
+ pandas
+ tqdm (completely obsolete but helpful while working with huge star files [shows progress])

## The best setup:

Install own python enviroment:
```
python3 -m venv new-env
```
Activate enviroment:
```
source new-env/bin/activate
```
Install packages:
```
pip install numpy matplotlib pandas gemmi tqdm
```
Run script from command line (convert to Relion 3.0 star):
```
python Relion_simple_star_parser.py --i run_data.star --convert
```
Plot two columns:
```
python Relion_simple_star_parser.py --i run_it025_data.star --plot _rlnMaxValueProbDistribution --plot _rlnNrOfSignificantSamples
```
![alt text](https://github.com/dzyla/Relion_star_parser/blob/master/histogram.png
)


Plot two columns, where 1 is the particle index:
```
python Relion_simple_star_parser.py --i run_it025_data.star --plot index --plot _rlnDefocusU --plot_type line
```
![alt text](https://github.com/dzyla/Relion_star_parser/blob/master/defocus.png
)

# Additional tools (not implemented via CLI yet)

It is possible to remove, modify and filter data in the star file. For now the only way to do it is changing the script. The data from star file is kept as the DataFrame, so all DataFrame options should work here too. 

At the bottom of the script are some examples:

#### 1. Divide the selected column by 10
```
particles_data_['_rlnMaxValueProbDistribution'] = particles_data_['_rlnMaxValueProbDistribution'].astype(float) / 10
```
#### 2. Remove one column completely
```
particles_data_ = particles_data_.drop(columns = ['_rlnOpticsGroup'])
```
#### 3. Filter data based on values in columns
```
particles_data_ = particles_data[particles_data['_rlnNrOfSignificantSamples'].astype(float) <= 10]
particles_data_ = particles_data_[particles_data['_rlnMaxValueProbDistribution'].astype(float) >= 0.5]
```
#### At the end save the the data as star:
```
save_star(particles_data_)
```


