# Written by Erik Madsen
# erik@madsense.net

import numpy as np
from astropy.io import fits
from astropy import units
import astropy.coordinates as coord
import astropy.time as atime
#from ch_util import ch_hdf5
from ch_util import andata
import h5py, os, argparse
from glob import glob

parser = argparse.ArgumentParser(description="Convert a directory of CHIME hdf5 files into a PSRFITS file.")
#parser.add_argument('-p', '--pulsar', metavar='NAME', dest='psr_name', default='pulsar', help="Name of pulsar observed. This is just used as FITS header information and defaults to \'pulsar\' if omitted.")
parser.add_argument('-r', '--ra', metavar='RA', dest='psr_ra', required=True, help="Pulsar right ascension in hh:mm:ss.sss form.")
parser.add_argument('-d', '--dec', metavar='DEC', dest='psr_dec', required=True, help="Pulsar declination in +dd:mm:ss.sss form.")
parser.add_argument('-f', '--firstfile', metavar='N', dest='first_file', type=int, default=0, help="First file in directory to read, counting from 0.")
parser.add_argument('-n', '--numfiles', metavar='N', dest='how_many_files', type=int, default=0, help="Number of files in directory to read (0 for all files).")
parser.add_argument('-s', '--npersubint', metavar='N', dest='n_per_subint', type=int, default=360, help="Number of samples per PSRFITS subintegration.")
#parser.add_argument('-h', '--infilesperoutfile', metavar='N', dest='in_per_out', type=int, default=0, help="Number of hdf5 files per FITS file. 0 uses all hdf5 files for one fits file. A nonzero value can result in npersubint changing.")
parser.add_argument('-o', '--outfile', metavar='FILE', dest='out_fname', help="Name of the fits file produced (defaults to [hdf5_dir].fits in working directory).")
parser.add_argument('-i', '--intermediate', action='store_true', dest='use_intermediate_file', help="If this flag is included, a file [hdf5_dir]_intmd.npy will be output after reading in the data or, if it already exists, will be read instead of the hdf5 files. This is mainly intended for debugging.")
parser.add_argument('-c', '--correctsampling', action='store_true', dest='correct_sampling', help="If this flag is included, the average spacing between timestamps will be used rather than the sampling time reported in the hdf5 headers.")
parser.add_argument('hdf5_dir', metavar='hdf5_directory', help="The directory containing the hdf5 files to be converted.")
args = parser.parse_args()

#psr_name = args.psr_name
psr_ra = args.psr_ra
psr_dec = args.psr_dec
first_file = args.first_file
how_many_files = args.how_many_files
n_per_subint = args.n_per_subint
if args.out_fname is not None:
    if args.out_fname[-5:] != ".fits":
        out_fname = args.out_fname + ".fits"
    else:
        out_fname = args.out_fname
else:
    out_fname = args.hdf5_dir + '.fits'
use_intermediate_file = args.use_intermediate_file
hdf5_dir = args.hdf5_dir

# NOTE for astropy 0.4, this should be changed from coord.ICRS to coord.SkyCoord
# (I think the rest of the syntax is the same)
psr_pos = coord.ICRS(ra=psr_ra, dec=psr_dec, unit=(units.hourangle, units.degree))
chime_longitude = -119.624

print "Reading CHIME data..."
intmd_file = os.path.basename(os.path.normpath(hdf5_dir)) + "_intmd.npy"
if os.path.exists(intmd_file) and use_intermediate_file:
    print "Loading intermediate file %s" % intmd_file
    file_contents = np.load(intmd_file)
    data_shape = tuple(file_contents[:2].astype(int))
    start_end = file_contents[2:4]
    tres, tres_apparent = file_contents[4:6]
    data = file_contents[6:].reshape(data_shape)
    del file_contents
else:
#    all_files = ch_hdf5.get_filepaths(hdf5_dir, first_file, how_many_files)
    all_files = glob(hdf5_dir + "/*.h5*")
    all_files.sort()
    with h5py.File(all_files[first_file], mode='r') as f:
        tres, = f.attrs['fpga.int_period']
        samples_per_hdf5, = f.attrs['acq.frames_per_file']
    if how_many_files == 0: how_many_files = len(all_files[first_file:])
#    last_file = first_file + how_many_files - 1
    data = np.zeros((samples_per_hdf5*how_many_files, 8, 1024), dtype='float32')
    tstamps = np.zeros(samples_per_hdf5*how_many_files)
    #data, tstamps = [], []
    for ii in range(how_many_files):
#    while first_file <= last_file:
        print "%d file(s) remaining..." % (how_many_files - ii)
        data_read = andata.Reader(all_files[first_file + ii])
        data_read.select_prod_autos()
        data_piece = data_read.read()
#        chan_means = data_piece.vis.real.T.mean(axis=0)
#        chan_means_div = chan_means.copy()
#        chan_means_div[chan_means == 0.] = 1.
#        data.append(((data_piece.vis.real.T - chan_means/chan_means_div).sum(axis=1)))
        lendiff = samples_per_hdf5 - data_piece.vis.real.T.shape[0]
        if lendiff > 0:
            if how_many_files - ii == 1:
                data = data[:-lendiff]
                tstamps = tstamps[:-lendiff]
                data[samples_per_hdf5*ii:] = data_piece.vis.real.T
                tstamps[samples_per_hdf5*ii:] = data_read.time
            else:
                print "File %d of %d appears to have the wrong number of samples." % (ii+1, how_many_files)
        else:
            data[samples_per_hdf5*ii:samples_per_hdf5*(ii+1)] = data_piece.vis.real.T
#        data.append(data_piece.vis.real.T)
            tstamps[samples_per_hdf5*ii:samples_per_hdf5*(ii+1)] = data_read.time
#        tstamps.append(data_read.time)
#        first_file += 1
#        data_piece, tstamps = ((data_piece[0] - chan_means)/chan_means_div).sum(axis=1), data_piece[1]
#    data = np.concatenate(data)
#    tstamps = np.concatenate(tstamps)

    chan_means = data.mean(axis=0)
    chan_means_div = chan_means.copy()
    chan_means_div[chan_means == 0.] = 1.
    #data = ((data - chan_means)/chan_means_div).sum(axis=1)
    data -= chan_means
    data /= chan_means_div
    data = data.sum(axis=1)

#    # In case last hdf5 file was shorter than the rest
#    tstamp_zero, = np.where(tstamps == 0)
#    tstamp_zero = tstamp_zero[tstamp_zero >= len(tstamps) - samples_per_hdf5]
#    if len(tstamp_zero):
#        tstamps = tstamps[:tstamp_zero[0]]
#        data = data[:tstamp_zero[0]]

    start_end = np.array([tstamps[0], tstamps[-1]])
    print "Done."
    # find missing samples and replace them with 'nan' for now
    print "Inserting missing samples..."
    nsamps_per_samp = (np.rint(np.diff(tstamps) / tres) + 0.01).astype(int)
    tres_apparent = np.diff(tstamps)[nsamps_per_samp == 1].mean()
    print "Apparent average timestep (s): %.8f" % tres_apparent
    insert_before = np.where(nsamps_per_samp > 1)[0] + 1
    size_of_gap = nsamps_per_samp[insert_before - 1] - 1
    all_inserts = np.ones((size_of_gap.sum(), data.shape[1]), dtype=bool)
    print "Inserting %d samples into %d gaps." % (all_inserts.shape[0], len(insert_before))
    where_inserts = []
    for ii in range(len(size_of_gap)):
        where_inserts += size_of_gap[ii] * [insert_before[ii]]
#    data = np.insert(data, where_inserts, -all_inserts.astype(float), axis=0)
    data = np.insert(data, where_inserts, all_inserts.astype(float)*np.nan, axis=0)
    if use_intermediate_file:
        print "Saving intermediate file %s" % intmd_file
        np.save(intmd_file, np.concatenate((np.array(data.shape).astype(float), start_end, np.array([tres, tres_apparent]), data.flatten())))
print "Done."

if args.correct_sampling:
    print "Using apparent sampling time of %.7f s rather than %.7f s." % (tres_apparent, tres)
    tres = tres_apparent

obs_start_time = atime.Time(start_end[0], format='unix', precision=0)
obs_end_time = atime.Time(start_end[1], format='unix', precision=0)
obs_duration = (obs_end_time - obs_start_time).sec

fits_data = fits.HDUList()
prim_hdr = fits.Header()
tbl_hdr = fits.Header()

current_time = atime.Time.now()
current_time.precision = 0

t_subint = n_per_subint * tres
n_subints = data.shape[0]/n_per_subint
n_chan = data.shape[1]

print "n_subints: %d" % n_subints
print "n_chan: %d" % n_chan

def lst(t, longitude):
    """
    t is an astropy.time.Time object in UTC
    longitude in decimal degrees: negative = W, positive = E
    returns lst in seconds

    a future version of astropy will have a sidereal_time function, yay
    """
    h = (t.mjd - int(t.mjd)) * 24.
    d = t.jd - 2451545.
    d0 = d - h/24.
    gst = 6.697374558 + 0.06570982441908*d0 + 1.00273790935*h + 0.000026*d/36525.
    return ((gst + longitude/15.) * 3600.) % 86400.

start_lst = lst(obs_start_time, chime_longitude)

prim_hdr.add_comment(
    "FITS (Flexible Image Transport System) format defined in Astronomy and")
prim_hdr.add_comment(
    "Astrophysics Supplement Series v44/p363, v44/p371, v73/p359, v73/p365.")
prim_hdr.add_comment(
    "Contact the NASA Science Office of Standards and Technology for the")
prim_hdr.add_comment(
    "FITS Definition document #100 and other FITS information.")
prim_hdr.add_comment("")
prim_hdr['hdrver'] = '5.3'
prim_hdr['fitstype'] = 'PSRFITS'
prim_hdr['date'] = current_time.isot
prim_hdr['observer'] = ''
prim_hdr['projid'] = ''
prim_hdr['telescop'] = 'CHIME'
#prim_hdr['ant_x'] = -2059164.942
#prim_hdr['ant_y'] = -3621108.403
#prim_hdr['ant_z'] = 4814432.276
prim_hdr['frontend'] = 'CHIME'
prim_hdr['nrcvr'] = 1
prim_hdr['fd_poln'] = 'LIN'
prim_hdr['fd_hand'] = 1
prim_hdr['fd_sang'] = 0.
prim_hdr['fd_xyph'] = 0.
prim_hdr['backend'] = 'CHIME'
prim_hdr['beconfig'] = 'N/A'
prim_hdr['be_phase'] = 0
prim_hdr['be_dcc'] = 0
prim_hdr['be_delay'] = 0.
prim_hdr['tcycle'] = 0.
prim_hdr['obs_mode'] = 'SEARCH'
prim_hdr['date-obs'] = obs_start_time.isot
prim_hdr['obsfreq'] = 600.1953125
prim_hdr['obsbw'] = -400.
prim_hdr['obsnchan'] = n_chan
prim_hdr['chan_dm'] = 0.
prim_hdr['src_name'] = os.path.basename(os.path.normpath(hdf5_dir))
prim_hdr['coord_md'] = 'J2000'
prim_hdr['equinox'] = 2000.
prim_hdr['ra'] = psr_ra
prim_hdr['dec'] = psr_dec
prim_hdr['bmaj'] = 0.
prim_hdr['bmin'] = 0.
prim_hdr['bpa'] = 0.
prim_hdr['stt_crd1'] = psr_ra
prim_hdr['stt_crd2'] = psr_dec
prim_hdr['trk_mode'] = 'TRACK'
prim_hdr['stp_crd1'] = psr_ra
prim_hdr['stp_crd2'] = psr_dec
prim_hdr['scanlen'] = np.round(obs_duration)
prim_hdr['fd_mode'] = 'FA'
prim_hdr['fa_req'] = 0.
prim_hdr['cal_mode'] = 'OFF'
prim_hdr['cal_freq'] = 0.
prim_hdr['cal_dcyc'] = 0.
prim_hdr['cal_phs'] = 0.
prim_hdr['stt_imjd'] = int(obs_start_time.mjd)
prim_hdr['stt_smjd'] = int(obs_start_time.mjd % 1 * 86400.)
prim_hdr['stt_offs'] = obs_start_time.mjd % 1 * 86400. % 1
prim_hdr['stt_lst'] = start_lst

prim_hdu = fits.PrimaryHDU(header=prim_hdr)

data = data[:n_subints*n_per_subint].reshape((n_subints, n_per_subint, n_chan))

print "Converting to 8-bit unsigned integers..."

#chan_means = np.ma.array(data, mask=(data < 0)).mean(axis=1)
chan_means = np.ma.array(data, mask=np.isnan(data)).mean(axis=1)
if chan_means.mask.any():
    print "Warning: at least one subint is entirely dummy samples. You might want to use longer subints."
#data[data < 0] = np.tile(chan_means.data, n_per_subint).reshape(data.shape)[data < 0]
data[np.isnan(data)] = np.tile(chan_means.data, n_per_subint).reshape(data.shape)[np.isnan(data)]
#chan_means_divide = chan_means.data.copy()
#chan_means_divide[chan_means.data == 0.] = 1.
#data = (data - np.tile(chan_means.data, n_per_subint).reshape(data.shape)) / np.tile(chan_means_divide, n_per_subint).reshape(data.shape)

data -= np.tile(np.min(data, axis=1), n_per_subint).reshape(data.shape)
max_data = np.max(data, axis=1)
max_data[max_data == 0.] = 1.
data *= np.tile(255./max_data, n_per_subint).reshape(data.shape)

tbl_hdr['pol_type'] = "AA"
tbl_hdr['tbin'] = tres
tbl_hdr['nchan'] = n_chan
tbl_hdr['chan_bw'] = -400./n_chan
tbl_hdr['nchnoffs'] = 0
tbl_hdr['npol'] = 1
tbl_hdr['nsblk'] = n_per_subint
tbl_hdr['nbits'] = 8

print "Defining FITS columns..."
weights = np.ones((n_subints, n_chan))
weights[:,992:] = 0. # this effectively removes the zero signal channels
offs_sub = np.arange(0, n_subints*t_subint, t_subint) + 0.5*t_subint
lst_sub = np.array([lst(t, chime_longitude) for t in obs_start_time + atime.TimeDelta(offs_sub, format='sec')])
the_columns = [
    fits.Column(name="TSUBINT", format='1D', unit='s', array=np.ones(n_subints)*t_subint),
    fits.Column(name="OFFS_SUB", format='1D', unit='s', array=offs_sub),
    fits.Column(name="LST_SUB", format='1D', unit='s', array=lst_sub),
    fits.Column(name="RA_SUB", format='1D', unit='deg', array=np.ones(n_subints)*psr_pos.ra.degree),
    fits.Column(name="DEC_SUB", format='1D', unit='deg', array=np.ones(n_subints)*psr_pos.dec.degree),
    fits.Column(name="GLON_SUB", format='1D', unit='deg', array=np.ones(n_subints)*psr_pos.galactic.l.degree),
    fits.Column(name="GLAT_SUB", format='1D', unit='deg', array=np.ones(n_subints)*psr_pos.galactic.b.degree),
    fits.Column(name="FD_ANG", format='1E', unit='deg', array=np.zeros(n_subints)),
    fits.Column(name="POS_ANG", format='1E', unit='deg', array=np.zeros(n_subints)),
    fits.Column(name="PAR_ANG", format='1E', unit='deg', array=np.zeros(n_subints)),
    fits.Column(name="TEL_AZ", format='1E', unit='deg', array=np.zeros(n_subints)),
    fits.Column(name="TEL_ZEN", format='1E', unit='deg', array=np.zeros(n_subints)),
    fits.Column(name="DAT_FREQ", format='%dE'%n_chan, unit='MHz', array=np.tile(np.linspace(800, 400, n_chan, endpoint=False) - 400./n_chan/2., n_subints).reshape(n_subints, n_chan)),
    fits.Column(name="DAT_WTS", format='%dE'%n_chan, array=weights),
    fits.Column(name="DAT_OFFS", format='%dE'%n_chan, array=np.zeros((n_subints, n_chan))),
    fits.Column(name="DAT_SCL", format='%dE'%n_chan, array=np.ones((n_subints, n_chan))),
    fits.Column(name="DATA", format=str(n_chan*n_per_subint) + 'B', dim='(%d,1,%d)' % (n_chan, n_per_subint), array=data.reshape((n_subints, n_per_subint*n_chan))),
]

print "Building table..."
table_hdu = fits.BinTableHDU(fits.FITS_rec.from_columns(the_columns), name="subint", header=tbl_hdr)

print "Writing fits file..."
fits_data.append(prim_hdu)
fits_data.append(table_hdu)
fits_data.writeto(out_fname, clobber=True)

print "Wrote file %s." % out_fname
