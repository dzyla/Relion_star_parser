import numpy as np
import pandas as pd
from gemmi import cif

'''
A simple STAR file parser based on the gemmi library.

Script written by Dawid Zyla, ETH Zurich.

For now script allows to plot selected columns and convert Relion 3.1 star files to 3.0. Outside the CLI there is possibility to
modify STAR file by removing columns, modifying column or filtering rows based on selected values. This however has to be done manually for now. 
At some point there will be some features added available from command line, such as:
* removing obsolete columns
* column modification by simple arithmetic equations
* row filtering based on the column values (e.g. _rlnNrOfSignificantSamples, _rlnMaxValueProbDistribution or just _rlnClassNumber)
* exporting star coordinates to cistem and cryosparc
* opening and modification of non data star files (selecting blocks, loops etc)
* saving Relion 3.1+ data star files (or star files with multiple blocks)
* maybe some more ideas
'''


def save_star(dataframe_, filename='out.star'):
    out_doc = cif.Document()
    out_particles = out_doc.add_new_block('particles', pos=-1)

    # Row number is required for the column names to save the STAR file e.g. _rlnNrOfSignificantSamples #33
    column_names = dataframe_.columns
    column_names_to_star = []
    for n, name in enumerate(column_names):
        column_names_to_star.append(name + ' #{}'.format(n + 1))

    loop = out_particles.init_loop('', column_names_to_star)
    data_rows = dataframe_.to_numpy().astype(str).tolist()

    for row in data_rows:
        loop.add_row(row)

    out_doc.write_file(filename)
    print('File "{}" saved.'.format(filename))


def save_star_31(dataframe_optics, dataframe_particles, filename='out.star'):
    # For now only Relion star 3.1+ can be saved as 3.1 star. Adding optics will be implemented later.

    out_doc = cif.Document()
    out_particles = out_doc.add_new_block('optics', pos=-1)

    # Row number is required for the column names to save the STAR file e.g. _rlnNrOfSignificantSamples #33
    column_names = dataframe_optics.columns
    column_names_to_star = []
    for n, name in enumerate(column_names):
        column_names_to_star.append(name + ' #{}'.format(n + 1))

    loop = out_particles.init_loop('', column_names_to_star)
    data_rows = dataframe_optics.to_numpy().astype(str).tolist()

    # save optics loop
    for row in data_rows:
        loop.add_row(row)

    out_particles = out_doc.add_new_block('particles', pos=-1)

    column_names = dataframe_particles.columns
    column_names_to_star = []
    for n, name in enumerate(column_names):
        column_names_to_star.append(name + ' #{}'.format(n + 1))

    loop = out_particles.init_loop('', column_names_to_star)
    data_rows = dataframe_particles.to_numpy().astype(str).tolist()

    # save particles loop
    for row in data_rows:
        loop.add_row(row)

    out_doc.write_file(filename)
    print('File "{}" saved.'.format(filename))


def convert_optics(optics_data_):
    # used for saving Relion 3.1 files with optics groups.
    # Changes the dict so values are list now.

    for key in optics_data_.keys():
        optics_data_[key] = [optics_data_[key]]

    return optics_data_


def convert_new_to_old(dataframe_, optics_group, filename, magnification='100000'):
    if optics_group == {}:
        print('File is already in Relion 3.0 format. No conversion needed!')
        quit()

    # change the Origin from Angstoms to pixels
    dataframe_['_rlnOriginXAngst'] = dataframe_['_rlnOriginXAngst'].astype(float) / optics_group[
        '_rlnImagePixelSize'].astype(float)
    dataframe_['_rlnOriginYAngst'] = dataframe_['_rlnOriginYAngst'].astype(float) / optics_group[
        '_rlnImagePixelSize'].astype(float)
    dataframe_ = dataframe_.rename(columns={'_rlnOriginXAngst': '_rlnOriginX', '_rlnOriginYAngst': '_rlnOriginY'})

    # add columns which are in the optics group
    dataframe_['_rlnVoltage'] = np.zeros(dataframe_.shape[0]) + optics_group['_rlnVoltage'].astype(float)
    dataframe_['_rlnSphericalAberration'] = np.zeros(dataframe_.shape[0]) + optics_group[
        '_rlnSphericalAberration'].astype(float)
    dataframe_['_rlnDetectorPixelSize'] = np.zeros(dataframe_.shape[0]) + optics_group['_rlnImagePixelSize'].astype(
        float)
    dataframe_['_rlnMagnification'] = np.zeros(dataframe_.shape[0]) + int(magnification)
    dataframe_['_rlnSphericalAberration'] = np.zeros(dataframe_.shape[0]) + optics_group[
        '_rlnSphericalAberration'].astype(float)

    # remove non used columns
    for tag in ['_rlnOpticsGroup', '_rlnHelicalTrackLengthAngst']:
        try:
            dataframe_ = dataframe_.drop(columns=[tag])
        except:
            pass

    # Row number is required for the column names
    column_names = dataframe_.columns
    column_names_to_star = []
    for n, name in enumerate(column_names):
        column_names_to_star.append(name + ' #{}'.format(n + 1))

    out_doc = cif.Document()
    out_particles = out_doc.add_new_block('', pos=-1)

    loop = out_particles.init_loop('', column_names_to_star)

    # to save cif all list values must be str
    data_rows = dataframe_.to_numpy().astype(str).tolist()

    for row in data_rows:
        loop.add_row(row)

    out_name = filename.replace('.star', '_v30.star')

    out_doc.write_file(out_name)
    print('File "{}" saved.'.format(out_name))


def parse_star(file_path):
    import tqdm

    doc = cif.read_file(file_path)

    optics_data = {}

    # 3.1 star files have two data blocks Optics and particles
    _new_star_ = True if len(doc) == 2 else False

    if _new_star_:
        print('Found Relion 3.1+ star file.')

        optics = doc[0]
        particles = doc[1]

        for item in optics:
            for optics_metadata in item.loop.tags:
                value = optics.find_loop(optics_metadata)
                optics_data[optics_metadata] = np.array(value)[0]

    else:
        print('Found Relion 3.0 star file.')
        particles = doc[0]

    particles_data = pd.DataFrame()

    print('Reading star file:')
    for item in particles:
        for particle_metadata in tqdm.tqdm(item.loop.tags):
            # If don't want to use tqdm uncomment bottom line and remove 'import tqdm'
            # for particle_metadata in item.loop.tags:
            loop = particles.find_loop(particle_metadata)
            particles_data[particle_metadata] = np.array(loop)

    return optics_data, particles_data


def parse_star_selected_columns(file_path, col1_name, col2_name):
    doc = cif.read_file(file_path)

    optics_data = {}

    # 3.1 star files have two data blocks Optics and particles
    _new_star_ = True if len(doc) == 2 else False

    if _new_star_:
        print('Found Relion 3.1+ star file.')

        optics = doc[0]
        particles = doc[1]

        for item in optics:
            for optics_metadata in item.loop.tags:
                value = optics.find_loop(optics_metadata)
                optics_data[optics_metadata] = np.array(value)[0]

    else:
        print('Found Relion 3.0 star file.')
        particles = doc[0]

    particles_data = pd.DataFrame()

    print('Reading star file:')

    for particle_metadata in [col1_name, col2_name]:
        loop = particles.find_loop(particle_metadata)
        particles_data[particle_metadata] = np.array(loop)

    return optics_data, particles_data


def plot_columns(particles_data, col1_name, col2_name, plot_type='hist'):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    print('Plotting {} on x and {} on y'.format(col1_name, col2_name))

    if col1_name != 'index' and col2_name != 'index':
        x_data = np.array(particles_data[col1_name].astype(float))
        y_data = np.array(particles_data[col2_name].astype(float))

    elif col1_name == 'index':
        x_data = np.arange(1, particles_data.shape[0] + 1, 1)
        y_data = np.array(particles_data[col2_name].astype(float))

    elif col2_name == 'index':
        y_data = np.arange(0, particles_data.shape[0], 1)
        x_data = np.array(particles_data[col1_name].astype(float))

    if plot_type == 'hist':
        plt.hist2d(x_data, y_data, cmap='Blues', bins=50, norm=LogNorm())
        clb = plt.colorbar()
        clb.set_label('Number of particles')

    elif plot_type == 'line':
        plt.plot(x_data, y_data)

    elif plot_type == 'scat':
        plt.scatter(x_data, y_data, cmap='Blues')

    plt.xlabel(col1_name)
    plt.ylabel(col2_name)

    plt.show()


'''----------------------------------------------------------------------------------------------------'''

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(
        description='Modify and plot Relion *.star files and convert Relion 3.1+ star to 3.0')
    parser.add_argument('--i', type=str, help='File path')
    parser.add_argument('--convert', action='store_true', default=False,
                        help='Convert Relion 3.1+ data.star to 3.0')
    parser.add_argument('--plot', action='append', dest='plot',
                        default=[],
                        help='plot two columns such as: --plot _rlnMaxValueProbDistribution --plot _rlnCtfFigureOfMerit. To plot the the value against index of the particle use --plot index --plot _rlsYOURMETA')
    parser.add_argument('--plot_type', dest='plot_type',
                        default='hist',
                        help='Plot type: histogram (hist) (default), scatter (scat) or XY plot (line) --plot_type scat or --plot_type hist --plot_type line')

    args = parser.parse_args()

    if not args.plot:
        try:
            path = args.i
            optics_data, particles_data = parse_star(path)
        except:
            # Here for the run without the command line
            quit()
            path = 'run_data.star'
            optics_data, particles_data = parse_star(path)

        print('Raw star file dimentions: {}'.format(particles_data.shape))

    else:
        try:
            path = args.i
            optics_data, particles_data = parse_star_selected_columns(path, args.plot[0], args.plot[1])
        except:
            quit()
            path = 'run_it025_data.star'
            optics_data, particles_data = parse_star_selected_columns(path, args.plot[0], args.plot[1])

    particles_data_ = particles_data

    if args.convert:
        convert_new_to_old(particles_data_, optics_data, args.i)

    if len(args.plot) == 2:
        try:
            plot_columns(particles_data, args.plot[0], args.plot[1], args.plot_type)

        except KeyError:
            print('Column not found. Check spelling or your star file')
            quit()

        except ValueError:
            print('Column is a text column. Not plottable.')


    elif len(args.plot) == 1 or len(args.plot) > 2:
        print('Only 2 columns are accepted.')
        quit()

    '''
    Modify star columns
    Uncomment the line if you want to change your star file, including:
    * changing the values in the column (add, subtract, multiply or divide)
    * select only particles (rows) which have specific values
    * Change column names
    * plot selected columns
    '''

    # Select rows with specific values
    # particles_data_ = particles_data[particles_data['_rlnNrOfSignificantSamples'].astype(float) <= 10]
    # particles_data_ = particles_data_[particles_data['_rlnMaxValueProbDistribution'].astype(float) >= 0.5]

    # Remove selected columns
    # particles_data_ = particles_data_.drop(columns = ['_rlnOpticsGroup'])

    # Modify values of the column
    # particles_data_['_rlnMaxValueProbDistribution'] = particles_data_['_rlnMaxValueProbDistribution'].astype(float) * 2

    # print('Modified star file dimentions: {}.\nRemoved {} particles'.format(particles_data_.shape, particles_data.shape[0]-particles_data_.shape[0]) )

    # Save modified file
    # save_star(particles_data_)
