# Relion_star_parser
A simple Relion .star file parser which allows Relion star file plotting, modification and conversion

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


