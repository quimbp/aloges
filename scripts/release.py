import json
import sys
import numpy
import datetime

# Sample over time, add a release point from each trajectory
# every SAMPLE_PERIOD hours
#
SAMPLE_OVER_TIME = False
SAMPLE_PERIOD = 24 # Hours, time between new releases.

# The release point can be retrieved from the point associated
# with the event 0 or the first point of the trajectory
# If RELEASE_INITIAL is False, then we read the 
# LineString data and get the first point that 
# after the RELEASE_DATE
#
# If RELEASE_DATE is set to None, but RELEASE_INITIAL
# is set to True, the the RELEASE_DATE is calcualted
# automatically.
#
RELEASE_INITIAL = True
RELEASE_IS_EVENT0 = True

RELEASE_DATE = None
RELEASE_DATE = datetime.datetime(2022,2,3,0,30,0) 

# Release Depth
#
RELEASE_DEPTH = 0.0
RELEASE_DEPTH_SVP = 15.0
RELEASE_DEPTH_DUMMYV = 0.5
RELEASE_DEPTH_SURFACE = 0.5
RELEASE_DEPTH_LIFERAFT = 0.0
 
# The output release filename:
#
RELEASE_FILENAME = 'release1.dat'

# Read argument options
# The first argument (i.e. filelist[0]) is the name of the script
#
filelist = list(sys.argv)
del filelist[0]

XP = []
YP = []
DP = []
TP = []

print('filelist: ', filelist)

for filename in filelist:

  print('\n %s' % filename)
  print(' --------------------------------')
  with open(filename,'r') as f:
    data = json.loads(f.read())

  FEATURES = data['features']

  for F in FEATURES:

    Ftype = F['geometry']['type']

    if Ftype == 'Point':

      if int(F['properties']['event']) == 0:
        particle_type = str(F['properties']['source'])
        experiment = str(F['properties']['exp'])
        x0 = F['geometry']['coordinates']
        d0 = F['properties']['time']['data'][0]
        d0 = d0[0:19]
        D0 = datetime.datetime.strptime(d0,'%Y-%m-%dT%H:%M:%S')
        D0N = datetime.datetime(D0.year,D0.month,D0.day,D0.hour,0,0)
        D0N += datetime.timedelta(hours=SAMPLE_PERIOD)
        print('  > Particle type:', particle_type)
        print('  > Experiment   :', experiment)
        print('  > Event 0 point         : %8.4f, %8.4f,  %s ' % (x0[0], x0[1], D0))

        # Append thes values to the Event 0 tables:
        #
        XE0 = x0[0]
        YE0 = x0[1]
        DE0 = D0

    else:

      # Here retrieve the first position and date of the Linstring:
      L_coords = F['geometry']['coordinates']
      L_time = F['properties']['time']['data']
      L_dates = []
      for tt in L_time:
        L_dates.append(datetime.datetime.strptime(tt[0:19],'%Y-%m-%dT%H:%M:%S'))
      
      x0 = L_coords[0]
      d0 = L_dates[0]
      print('  > First Linestring point: %8.4f, %8.4f,  %s ' % (x0[0], x0[1], d0))
    
      XL0 = x0[0]
      YL0 = x0[1]
      DL0 = D0

      x1 = L_coords[-1]
      d1 = L_dates[-1]
      D1P = datetime.datetime(d1.year,d1.month,d1.day,d1.hour,0,0)
      #D1P -= datetime.timedelta(hours=SAMPLE_PERIOD)
      print('  > Last  Linestring point: %8.4f, %8.4f,  %s ' % (x1[0], x1[1], d1))

      TP.append(particle_type)
      if RELEASE_INITIAL:
        I0 = 0
        if RELEASE_IS_EVENT0:
          print('     Adding event 0')
          XP.append(XE0)
          YP.append(YE0)
          DP.append(DE0)
        else:
          print('     Adding linestring first value')
          XP.append(XL0)
          YP.append(YL0)
          DP.append(DL0)
      else:
        I0 = -1
        for i in range(len(L_dates)):
          if L_dates[i] >= RELEASE_DATE:
            I0 = i
            XO0 = L_coords[i][0]
            YO0 = L_coords[i][1]
            DO0 = L_dates[i]
            break
        if I0 < 0:
          print('\nDate: ', RELEASE_DATE,' not found') 
          quit()
        else:
          print('     Adding linestring value at',I0,' of date', DO0,' :: ',XO0, YO0)
          XP.append(XO0)
          YP.append(YO0)
          DP.append(DO0)
          D0N = RELEASE_DATE
            
      


      if SAMPLE_OVER_TIME:

        #nperiods = int((D1P - D0N).total_seconds() / (SAMPLE_PERIOD*3600.0))
        nperiods = int((d1 - D0N).total_seconds() / (SAMPLE_PERIOD*3600.0))
        print('  > Sampling period [h]   : %8.4f' % (SAMPLE_PERIOD))
        print('  > Number of periods     : ', nperiods)

        DI = D0N
        DI += datetime.timedelta(hours=SAMPLE_PERIOD)
        for h in range(nperiods+1):

          for i in range(I0,len(L_dates)):
            if L_dates[i] >= DI:
              print('     Adding linestring value at',i,' of date', L_dates[i], ' :: ', L_coords[i][0],L_coords[i][1])
              XP.append(L_coords[i][0])
              YP.append(L_coords[i][1])
              DP.append(L_dates[i])
              break

          DI += datetime.timedelta(hours=SAMPLE_PERIOD)
      

if RELEASE_DATE is None:
  TMIN = numpy.min(DP)
  yy = TMIN.year
  mm = TMIN.month
  dd = TMIN.day
  release_date = datetime.datetime(yy,mm,dd,0,0,0)
else:
  release_date = RELEASE_DATE
  
string = datetime.datetime.strftime(release_date,'%Y-%m-%dT%H:%M:%S')

print('\nRelease date:', release_date)

# Remove duplicated positions:
#
Zall = []
for i in range(len(XP)):
  Zall.append([XP[i],YP[i],DP[i],TP[i]])

Z = []
for z in Zall:
  if z not in Z:
    Z.append(z)


with open(RELEASE_FILENAME,'w') as f:
  f.write('# Release date: '+string+'\n')
  f.write('# lon       lat       depth     secs\n')
  for z in Z:
    if z[3] == 'liferaft':
      DEPTH = RELEASE_DEPTH_LIFERAFT
    elif z[3] == 'dr_surface':
      DEPTH = RELEASE_DEPTH_SURFACE
    elif z[3] == 'dummy_v':
      DEPTH = RELEASE_DEPTH_DUMMYV
    elif z[3] == 'dr_svp':
      DEPTH = RELEASE_DEPTH_SVP
    else:
      DEPTH = RELEASE_DEPTH
    f.write('%9.6f   %9.6f   %9.3f  %7.0f\n' % (z[0], z[1], DEPTH, (z[2]-release_date).total_seconds()))
