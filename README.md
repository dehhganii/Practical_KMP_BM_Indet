# Practical_KMP_BM_Indet

### How to use the code 
You can create files with using
`make all` command.
<br><br><br>
<br><br>
### Generate new strings
For generating new strings you can use the command <br>
`python3 create_strings.py {n} {k1} {m} {k2} {sigma}`
Which are:

  * n = text size
  * k1 = number of indeterminate letters in text
  * m = pattern size
  * k2 = number of indeterminate letters in text
  * sigma = alphabet size
<br><br><br>
<br><br>
### For using the generated file you should run these commands
* Brute force  `./brute_force {alphabet_size} {text_file_name} {pattern_file_name}`
* Kmp indet  `./kmp_indet {alphabet_size} {text_file_name} {pattern_file_name}`
* Boyer moore indet  `./bm_indet {alphabet_size} {text_file_name} {pattern_file_name}`
<br><br><br>
<br><br>
##### In the example text and pattern in the repository the alphabet size is **4**

#### If you want to run all three algorithms at the same time you can use `run_all_algorithms` like this:
 `./run_all_algorithms {alphabet_size}`
 
*just remember to do `chmod +x run_all_algorithms` to give execution permission to the file*
