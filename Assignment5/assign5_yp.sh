#Create parent directory
mkdir Assign5

#switch to new directory
cd Assign5

#Create 500 subdirectories which all include a file with 5 lines in each file
for num in {001..500};
do
mkdir dir$num
cd dir$num
printf '%s\n' 'line1' 'line2' 'line3' 'line4' 'line5' > file.txt	
cd ..
done
