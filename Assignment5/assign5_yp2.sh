#switch to parent directory
cd /auto/pmd-01/ypatel/Assign5

#If the subdirectory ends in the number 1, the fourth line of text should read: "eat beets"
for dir in -d -- *[1]/; 
do
cd "$dir"
awk 'NR==4 {$0="eat beets"} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done

#If the subdirectory ends in the number 4, the fourth line of text should read: "squash are great"
for dir in -d -- *[4]/; 
do
cd "$dir"
awk 'NR==4 {$0="squash are great"} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done

#If the subdirectory ends in the number 5, the fourth line of text should read: "dogs are better than cats"
for dir in -d -- *[5]/; 
do
cd "$dir"
awk 'NR==4 {$0="dogs are better than cats"} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done

#If the subdirectory ends in the number 7, the fourth line of text should read: "hello world"
for dir in -d -- *[7]/; 
do
cd "$dir"
awk 'NR==4 {$0="hello world"} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done

#If the subdirectory ends in the number 0, the fourth line of text should be a phrase of your choosing.
for dir in -d -- *[0]/; 
do
cd "$dir"
awk 'NR==4 {$0="eureka i did it"} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done

#If the subdirectory ends in the number 2, 3, 6, 8, or 9, the fourth line of text should be blank. 
for dir in -d -- *[2,3,6,8,9]/; 
do
cd "$dir"
awk 'NR==4 {$0=" "} 1' file.txt > file_tmp.txt && mv file_tmp.txt file.txt
cd /auto/pmd-01/ypatel/Assign5
done
