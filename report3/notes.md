# Part 2 Notes:

navigate to `data` file
```python
# Copy all directories
cp ../bias/b
cp ../bias/bias.sdf
cp ../flat/balance.sdf
cp ../flat.zero.sdf
cp ../flat/unit.sdf
ls
pplot rXXX.sdf
gai
gaia r
gaia rXXX.sdf
ls
./debias
# chech it works by replotting
pplot rXXX.sdf
```
Graph should be from 300 to 4000 roughly with spikes throughout the spectrum

```python
emacs noise &
```
Set image to desired file. This is going to perform a subtraction of all the spectra.

1. Choose image (spectrum of one of our stars)
```python     
gaia r 
gaia rXXX.sdf
```
- If the name is something like SP09043-441 it's good?
Search: detector ==> should be: EEV12

- Search WHT EEV12 on Google and get the Gain value from the first link. This is in (e/ADU). Put this in the photon = x value in the noise window

- Readout Speed: SLOW (from Gaia), Also on the internet. Assign this to readout = xx

- ystart and yend are fixed: 140 - 

```python
./n
./noise
# do you want frame plot? [y]
y
# Plot fit? [y]
y
y
ent
ent
ent
ent
ent
c # to continue
```
on left side of line: 
```
c
c
```
on right side of line
```
c
c
```
`q`

Look at average RMS. Should be close to one. 3-4 is wrong

```
y 
y
y 
y
```

Look at reduced Chi-squared of fit. Should be near 1

Shift is -0.42. Should be very small
```
ls
emacs datafile &
ls r* > reduce.lis
emacs reduce.lis
```
Change xstart, xend, ystart, yend, trackstart, track_end, skymov_slo

```
reduce dat
reduce -t datafile
ls zzz_datatypes
emacs zzz_datatypes &
```
Data means spectrum of a star, ARC means something else

```
reduce datafile
```
Now it should finish
```python
ls
# Should be tons of files

pplot rxxxopt.sdf # Optimal Structure
# Should give a crazy graph

pplot rxxxnor.sdf # Normal Region

pplot rxxxskyo.sdf # sky Optimal Structure

pplot rxxxskyn.sdf # Sky Normal Region

pplot rxxxarco.sdf # ARC Optimal

# Each object should have one associated spectrum


```