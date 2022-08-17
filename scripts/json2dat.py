import json
import sys
import datetime
import numpy as np
from scipy import interpolate

INTERPOLATION_METHODS = {0: 'linear', 
                         1: 'nearest', 
                         2: 'nearest-up', 
                         3: ' zero',
                         4: 'slinear', 
                         5: 'quadratic',
                         6: 'cubic', 
                         7: 'previous', 
                         8: 'next' }

VERSION = "1.0"
NOW = datetime.datetime.now()

INTERPOL = False
INTERPOL_DT   = 30 * 60    # seconds
INTERPOL_KIND = INTERPOLATION_METHODS[0]

# Release Depth
#
RELEASE_DEPTH = 0.0
RELEASE_DEPTH_SVP = 15.0
RELEASE_DEPTH_DUMMYV = 0.5
RELEASE_DEPTH_SURFACE = 0.5
RELEASE_DEPTH_LIFERAFT = 0.0

args = list(sys.argv)
nargs = len(args)

WithHelp = False
for ff in args:
  if '--help' in ff:
    WithHelp = True

if nargs < 3 or WithHelp:
  print('\njson2dat: python json2dat [FROM=INITIAL_DATE] [DT=SECONDS] [KIND=INTERPOLATION_KIND] in out')
  print('     Converts a GEOJSON trajectory file to a tabular text file.')
  print('     Version 1.0, August 2022.')
  print('')
  print('     Options:')
  print('       INITIAL_DATE expressed as ISO format (E.g.: D0=2022-12-31T23:00:00).')
  print('       DT           is the interpolation time step in seconds (E.g. DT=1800).')
  print('       KIND         is the interpolation kind. Must be an integer from 0 to 8 (E.g. KIND=0).\n')
  print('     Interpolation options, INTERPOLATION_KIND:')
  print('       0            Linear interpolation.')
  print('       1            Nearest.')
  print('       2            Nearest-up.')
  print('       3            Zero-order splines.')
  print('       4            Linear splines.')
  print('       5            Quadratic splines.')
  print('       6            Cubic splines.')
  print('       7            Previous value.\n')
  print('     Example:\n')
  print('       python json.dat -interpol from=2022-03-15T03:00:00 dt=30 kind=0 filein.geojson fileout.dat\n')

  quit()

Ifile = args[-2]
Ofile = args[-1]

WithD0 = False
WithDT = False
WithKn = False
for ff in args:
  FF = ff.upper()

  if 'INTERPOL' in FF:
    INTERPOL = True
  if 'FROM' in FF:
    WithD0 = True
    Dini = ff
  if 'DT' in FF:
    WithDT = True
    DT = ff
  if 'KIND' in FF:
    WithKn = True
    kinda = ff

if WithD0:
  i = Dini.find('=') + 1
  INTERPOL_D0 = Dini[i:]
  INTERPOL_D0 = INTERPOL_D0[0:19]

if WithDT:
  i = DT.find('=') + 1
  INTERPOL_DT = float(int(DT[i:]))

if WithKn:
  i = kinda.find('=') + 1
  Kind = int(kinda[i:])
  INTERPOL_KIND = INTERPOLATION_METHODS[Kind]

# .........................................................
# ... Read the data:
# ...
print('\nReading file ', Ifile)
with open(Ifile,'r') as f:
  data = json.loads(f.read())


XP = []
YP = []
ZP = []
TP = []
DP = []

FEATURES = data['features']

for F in FEATURES:

  Ftype = F['geometry']['type']

  if Ftype == 'Point':

    if int(F['properties']['event']) == 0:
      try:
        particle_type = str(F['properties']['source'])
      except:
        particle_type = ""
      try:
        experiment = str(F['properties']['exp'])
      except:
        experiment = ""
      try:
        code_sn = str(F['properties']['code_sn'])
      except:
        code_sn = ""
      try:
        ptype = str(F['properties']['tipo'])
      except:
        ptype = ""
      x0 = F['geometry']['coordinates']
      d0 = F['properties']['time']['data'][0]
      D0 = datetime.datetime.strptime(d0[0:19],'%Y-%m-%dT%H:%M:%S')
      T0 = D0.timestamp()
      print('  > Particle type:', particle_type)
      print('  > Experiment   :', experiment)
      print('  > Event 0 point         : %8.4f, %8.4f,  %s ' % (x0[0], x0[1], D0))

      # Append thes values to the Event 0 tables:
      #
      XP.append(float(x0[0]))
      YP.append(float(x0[1]))
      DP.append(D0)
      TP.append(T0)

      if particle_type == 'liferaft':
        DEPTH = RELEASE_DEPTH_LIFERAFT
      elif particle_type == 'dr_surface':
        DEPTH = RELEASE_DEPTH_SURFACE
      elif particle_type == 'dummy_v':
        DEPTH = RELEASE_DEPTH_DUMMYV
      elif particle_type == 'dr_svp':
        DEPTH = RELEASE_DEPTH_SVP
      else:
        DEPTH = RELEASE_DEPTH

    else:
      x1 = F['geometry']['coordinates']
      d1 = F['properties']['time']['data'][0]
      D1 = datetime.datetime.strptime(d1[0:19],'%Y-%m-%dT%H:%M:%S')
      T1 = D1.timestamp()

  else:

    # Here retrieve the first position and date of the Linstring:
    L_coords = F['geometry']['coordinates']
    L_time = F['properties']['time']['data']
    L_dates = []
    for tt in L_time:
      L_dates.append(datetime.datetime.strptime(tt[0:19],'%Y-%m-%dT%H:%M:%S'))

    for i in range(len(L_time)):
      xy = L_coords[i]
      XP.append(float(L_coords[i][0]))
      YP.append(float(L_coords[i][1]))
      DP.append(L_dates[i])
      TP.append(np.float(L_dates[i].timestamp()))

if TP[0] == TP[1]:
  print('WARNING: The first two values coincide. Removing first record !')
  XP = XP[1:]
  YP = YP[1:]
  TP = TP[1:]
  DP = DP[1:]


XP = np.array(XP)
YP = np.array(YP)
TP = np.array(TP)
DP = np.array(DP)
# ...
# .........................................................


if ( not INTERPOL):
  with open(Ofile,'w') as F:
    F.write("# ALOGES JSON2DAT.PY\n")
    F.write("# -----------------------------------------------------------\n")
    F.write("# Version: %s\n" % (VERSION))
    F.write("# Runtime: %s\n" % (NOW.replace(microsecond=0).isoformat()))
    F.write("# Trajectory from: %s\n" % (Ifile))
    F.write("# Particle type  : %s\n" % (particle_type))
    F.write("# Experiment     : %s\n" % (experiment))
    F.write("# code_sn        : %s\n" % (code_sn))
    F.write("# type           : %s\n" % (ptype))
    F.write("# Event 0 point  : %8.4f, %8.4f,  %s \n" % (x0[0], x0[1], D0.replace(microsecond=0).isoformat()))
    F.write("# Event 1 point  : %8.4f, %8.4f,  %s \n" % (x1[0], x1[1], D1.replace(microsecond=0).isoformat()))
    F.write("# -----------------------------------------------------------\n")
    for i in range(len(DP)):
      F.write("%12.6f  %12.6f   %9.3f    %s\n" % (XP[i],YP[i], DEPTH, DP[i].replace(microsecond=0).isoformat()))
  print('Data saved in ', Ofile)
  quit()



# ----------------------------------------------------------------------
#                             INTERPOL
# ----------------------------------------------------------------------
# If we are here is because we need to interpolate the trajectory data:
#

dt = DP[2] - DP[1]
dt = dt.total_seconds()/60.0
dt = np.rint(dt/10.0)*10.0

if not WithD0:
  d0 = DP[0]
  if d0.minute < 30:
    DI0 = d0.replace(minute=0,second=0)
  elif d0.minute > 30:
    DI0 = d0.replace(hour=d0.hour+1,minute=0,second=0)
else:
  DI0 = datetime.datetime.strptime(INTERPOL_D0,'%Y-%m-%dT%H:%M:%S')
   

print(" ")
print('Interpolation D0    : ', DI0)
print('Interpolation dt    : ', INTERPOL_DT)
print('Interpolation method: ', INTERPOL_KIND)


TI0 = DI0.timestamp()
nsteps = int(np.rint((T1 - TI0)/INTERPOL_DT)+1)

TIME = np.array([TI0+i*INTERPOL_DT for i in range(nsteps)])
DATE = []
for tt in TIME:
  DATE.append(datetime.datetime.fromtimestamp(tt))

f = interpolate.interp1d(TP,XP,kind=INTERPOL_KIND,bounds_error=False, fill_value=-999.0)
XI = f(TIME)

f = interpolate.interp1d(TP,YP,kind=INTERPOL_KIND,bounds_error=False, fill_value=-999.0)
YI = f(TIME)


with open(Ofile,'w') as F:
  F.write("# ALOGES JSON2DAT.PY\n")
  F.write("# -----------------------------------------------------------\n")
  F.write("# Version: %s\n" % (VERSION))
  F.write("# Runtime: %s\n" % (NOW.replace(microsecond=0).isoformat()))
  F.write("# Trajectory from: %s\n" % (Ifile))
  F.write("# Particle type  : %s\n" % (particle_type))
  F.write("# Experiment     : %s\n" % (experiment))
  F.write("# code_sn        : %s\n" % (code_sn))
  F.write("# type           : %s\n" % (ptype))
  F.write("# Event 0 point  : %8.4f, %8.4f,  %s \n" % (x0[0], x0[1], D0.replace(microsecond=0).isoformat()))
  F.write("# Event 1 point  : %8.4f, %8.4f,  %s \n" % (x1[0], x1[1], D1.replace(microsecond=0).isoformat()))
  F.write("#           >> INTERPOLATED TRAJECTORY : %s <<\n" % (INTERPOL_KIND))
  F.write("# -----------------------------------------------------------\n")
  for i in range(len(TIME)):
    F.write("%12.6f  %12.6f   %9.3f    %s\n" % (XI[i],YI[i], DEPTH, DATE[i].replace(microsecond=0).isoformat()))

print('Interpolated output saved in ', Ofile)
