#
# This is a starter Automake file for ridgerunner created 3/22/7
#

AM_CFLAGS = -g -finline-functions -Wall

bin_PROGRAMS = ridgerunner
ridgerunner_SOURCES = ./src/display.c \
	./src/errors.c \
	./src/linklib_additions.c \
	./src/settings.c \
	./src/dlen.c \
	./src/ridgerunner_main.c \
        ./src/stepper.c \
	./src/dlen.h \
	./src/eqedge.h \
	./src/errors.h \
        ./src/settings.h \
	./src/stepper.h 

ridgerunner_LDADD = -lplCurve -ltsnnls 

# We include a copy of the knotplot distribution of loose knots for testing.

data_DIST = 0.1.vect 0.2.1.vect 0.3.1.vect 0.4.1.vect 10.1.vect  \
	10.10.vect 10.100.vect 10.101.vect 10.102.vect 10.103.vect  \
	10.104.vect 10.105.vect 10.106.vect 10.107.vect 10.108.vect  \
	10.109.vect 10.11.vect 10.110.vect 10.111.vect 10.112.vect  \
	10.113.vect 10.114.vect 10.115.vect 10.116.vect 10.117.vect  \
	10.118.vect 10.119.vect 10.12.vect 10.120.vect 10.121.vect  \
	10.122.vect 10.123.vect 10.124.vect 10.125.vect 10.126.vect  \
	10.127.vect 10.128.vect 10.129.vect 10.13.vect 10.130.vect  \
	10.131.vect 10.132.vect 10.133.vect 10.134.vect 10.135.vect  \
	10.136.vect 10.137.vect 10.138.vect 10.139.vect 10.14.vect  \
	10.140.vect 10.141.vect 10.142.vect 10.143.vect 10.144.vect  \
	10.145.vect 10.146.vect 10.147.vect 10.148.vect 10.149.vect  \
	10.15.vect 10.150.vect 10.151.vect 10.152.vect 10.153.vect  \
	10.154.vect 10.155.vect 10.156.vect 10.157.vect 10.158.vect  \
	10.159.vect 10.16.vect 10.160.vect 10.161.vect 10.162.vect  \
	10.163.vect 10.164.vect 10.165.vect 10.166.vect 10.17.vect  \
	10.18.vect 10.19.vect 10.2.vect 10.20.vect 10.21.vect  \
	10.22.vect 10.23.vect 10.24.vect 10.25.vect 10.26.vect  \
	10.27.vect 10.28.vect 10.29.vect 10.3.vect 10.30.vect  \
	10.31.vect 10.32.vect 10.33.vect 10.34.vect 10.35.vect  \
	10.36.vect 10.37.vect 10.38.vect 10.39.vect 10.4.vect  \
	10.40.vect 10.41.vect 10.42.vect 10.43.vect 10.44.vect  \
	10.45.vect 10.46.vect 10.47.vect 10.48.vect 10.49.vect  \
	10.5.vect 10.50.vect 10.51.vect 10.52.vect 10.53.vect  \
	10.54.vect 10.55.vect 10.56.vect 10.57.vect 10.58.vect  \
	10.59.vect 10.6.vect 10.60.vect 10.61.vect 10.62.vect  \
	10.63.vect 10.64.vect 10.65.vect 10.66.vect 10.67.vect  \
	10.68.vect 10.69.vect 10.7.vect 10.70.vect 10.71.vect  \
	10.72.vect 10.73.vect 10.74.vect 10.75.vect 10.76.vect  \
	10.77.vect 10.78.vect 10.79.vect 10.8.vect 10.80.vect  \
	10.81.vect 10.82.vect 10.83.vect 10.84.vect 10.85.vect  \
	10.86.vect 10.87.vect 10.88.vect 10.89.vect 10.9.vect  \
	10.90.vect 10.91.vect 10.92.vect 10.93.vect 10.94.vect  \
	10.95.vect 10.96.vect 10.97.vect 10.98.vect 10.99.vect  \
	2.2.1.vect 3.1.vect 4.1.vect 4.2.1.vect 5.1.vect  \
	5.2.1.vect 5.2.vect 6.1.vect 6.2.1.vect 6.2.2.vect  \
	6.2.3.vect 6.2.vect 6.3.1.vect 6.3.2.vect 6.3.3.vect  \
	6.3.vect 7.1.vect 7.2.1.vect 7.2.2.vect 7.2.3.vect  \
	7.2.4.vect 7.2.5.vect 7.2.6.vect 7.2.7.vect 7.2.8.vect  \
	7.2.vect 7.3.1.vect 7.3.vect 7.4.vect 7.5.vect  \
	7.6.vect 7.7.vect 8.1.vect 8.10.vect 8.11.vect  \
	8.12.vect 8.13.vect 8.14.vect 8.15.vect 8.16.vect  \
	8.17.vect 8.18.vect 8.19.vect 8.2.1.vect 8.2.10.vect  \
	8.2.11.vect 8.2.12.vect 8.2.13.vect 8.2.14.vect 8.2.15.vect  \
	8.2.16.vect 8.2.2.vect 8.2.3.vect 8.2.4.vect 8.2.5.vect  \
	8.2.6.vect 8.2.7.vect 8.2.8.vect 8.2.9.vect 8.2.vect  \
	8.20.vect 8.21.vect 8.3.1.vect 8.3.10.vect 8.3.2.vect  \
	8.3.3.vect 8.3.4.vect 8.3.5.vect 8.3.6.vect 8.3.7.vect  \
	8.3.8.vect 8.3.9.vect 8.3.vect 8.4.1.vect 8.4.2.vect  \
	8.4.3.vect 8.4.vect 8.5.vect 8.6.vect 8.7.vect  \
	8.8.vect 8.9.vect 9.1.vect 9.10.vect 9.11.vect  \
	9.12.vect 9.13.vect 9.14.vect 9.15.vect 9.16.vect  \
	9.17.vect 9.18.vect 9.19.vect 9.2.1.vect 9.2.10.vect  \
	9.2.11.vect 9.2.12.vect 9.2.13.vect 9.2.14.vect 9.2.15.vect  \
	9.2.16.vect 9.2.17.vect 9.2.18.vect 9.2.19.vect 9.2.2.vect  \
	9.2.20.vect 9.2.21.vect 9.2.22.vect 9.2.23.vect 9.2.24.vect  \
	9.2.25.vect 9.2.26.vect 9.2.27.vect 9.2.28.vect 9.2.29.vect  \
	9.2.3.vect 9.2.30.vect 9.2.31.vect 9.2.32.vect 9.2.33.vect  \
	9.2.34.vect 9.2.35.vect 9.2.36.vect 9.2.37.vect 9.2.38.vect  \
	9.2.39.vect 9.2.4.vect 9.2.40.vect 9.2.41.vect 9.2.42.vect  \
	9.2.43.vect 9.2.44.vect 9.2.45.vect 9.2.46.vect 9.2.47.vect  \
	9.2.48.vect 9.2.49.vect 9.2.5.vect 9.2.50.vect 9.2.51.vect  \
	9.2.52.vect 9.2.53.vect 9.2.54.vect 9.2.55.vect 9.2.56.vect  \
	9.2.57.vect 9.2.58.vect 9.2.59.vect 9.2.6.vect 9.2.60.vect  \
	9.2.61.vect 9.2.7.vect 9.2.8.vect 9.2.9.vect 9.2.vect  \
	9.20.vect 9.21.vect 9.22.vect 9.23.vect 9.24.vect  \
	9.25.vect 9.26.vect 9.27.vect 9.28.vect 9.29.vect  \
	9.3.1.vect 9.3.10.vect 9.3.11.vect 9.3.12.vect 9.3.13.vect  \
	9.3.14.vect 9.3.15.vect 9.3.16.vect 9.3.17.vect 9.3.18.vect  \
	9.3.19.vect 9.3.2.vect 9.3.20.vect 9.3.21.vect 9.3.3.vect  \
	9.3.4.vect 9.3.5.vect 9.3.6.vect 9.3.7.vect 9.3.8.vect  \
	9.3.9.vect 9.3.vect 9.30.vect 9.31.vect 9.32.vect  \
	9.33.vect 9.34.vect 9.35.vect 9.36.vect 9.37.vect  \
	9.38.vect 9.39.vect 9.4.1.vect 9.4.vect 9.40.vect  \
	9.41.vect 9.42.vect 9.43.vect 9.44.vect 9.45.vect  \
	9.46.vect 9.47.vect 9.48.vect 9.49.vect 9.5.vect  

