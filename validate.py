#
# This script works only with Python 3
#
import os
import re

#import sys
#import time
#import random



import argparse
from zipfile import ZipFile 


import create_initial_conditions
import Test


#import AnswerDatabase
#import __main__


#MaxParticlesInSequentialUpscalingStudies = 2000
MaxParticlesInSequentialUpscalingStudies = 200
#MaxParticlesInSequentialUpscalingStudies = 200


def step1(zip):
  f = zip.open("step-1.cpp")
  content = f.read()
  f = open('step-1.cpp', 'wb')
  f.write(content)
  f.close()
  print( "========" )
  print( " Step 1" )
  print( "========" )
  print( "Extracted step-1.cpp ... ok" )
       
  test = Test.Test( "step-1.cpp" )
  
  arguments = "icpc --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )

  arguments = create_initial_conditions.create_no_noise_grid_setup( 0, 1, 0.0001, 0.1, 0.1, [2,1,1] )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )

  search_pattern = "([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+), *([-+]?\d*\.\d+|\d+)"
  result = test.search_for_pattern_in_output(search_pattern)
  if result!="":
    print( "Got " + result + " as last output which I interprete to be a particle position ... ok (though that does not mean that the data is correct; that's something I don't validate here)" )
  else:
    print( "Last line of output should still be the one I used in the template ... failed" )


def step2(zip):
  f = zip.open("step-2.cpp")
  content = f.read()
  f = open('step-2.cpp', 'wb')
  f.write(content)
  f.close()
  print( "========" )
  print( " Step 2" )
  print( "========" )
  print( "Extracted step-2.cpp ... ok" )
       
  test = Test.Test( "step-2.cpp" )
  
  arguments = "icpc --std=c++0x -no-vec"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )

  particles_counts = [11,11,11]
  arguments = create_initial_conditions.create_random_grid_setup( 0, 10, 0.0001, 0.1, 0.01, particles_counts )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
  no_vec_time = test.runtime
  
  arguments = "icpc -O3 -xhost -fopenmp --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )

  particles_counts = [11,11,11]
  arguments = create_initial_conditions.create_random_grid_setup( 0, 10, 0.0001, 0.1, 0.01, particles_counts )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok" )
  else:
    print( "Run failed: " + test.last_output )
  with_vec_time = test.runtime
  
  if no_vec_time<=with_vec_time:
    print( "Code is slower with vectorisation, so you might want to tune it ... check" )
  else:
    print( "Code is already faster by a factor of " + str(no_vec_time/with_vec_time) + " through vectorisation but you might want to tune it further ... ok" )


def step3(zip):
  f = zip.open("step-3.cpp")
  content = f.read()
  f = open('step-3.cpp', 'wb')
  f.write(content)
  f.close()
  print( "========" )
  print( " Step 3" )
  print( "========" )
  print( "Extracted step-3.cpp ... ok" )
       
  test = Test.Test( "step-3.cpp" )
  
  arguments = "icpc -O3 -xhost -fopenmp --std=c++0x"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )

  arguments = create_initial_conditions.create_no_noise_grid_setup( 0, 1, 0.0001, 0.1, 0.1, [2,1,1] )
  print( "Run code with " + arguments )
  result = test.run( arguments )
  if result==1:
    print( "Run terminated successfully ... ok (but this does not mean that the outcome is correct)" )
  else:
    print( "Run failed: " + test.last_output )


def step4(zip):
  f = zip.open("step-4.cpp")
  content = f.read()
  f = open('step-4.cpp', 'wb')
  f.write(content)
  f.close()
  print( "========" )
  print( " Step 4" )
  print( "========" )
  print( "Extracted step-4.cpp ... ok" )
       
  test = Test.Test( "step-4.cpp" )
  
  arguments = "icpc -O3 -xhost --std=c++0x -fopenmp"
  result = test.compile( arguments )
  if result==1:
    print( "Compiled source code with " + arguments + " ... ok" )
  else:
    print( "Compiled source code with " + arguments + " ... failed: " + test.last_output )

  particles_counts = [11,11,11]
  arguments = create_initial_conditions.create_random_grid_setup( 0, 10, 0.0001, 0.1, 0.01, particles_counts )
  print( "Run code with " + str(particles_counts[0]*particles_counts[1]*particles_counts[2]) + " particles" )

  print( "Test for one core" )
  result = test.run( arguments, {"OMP_NUM_THREADS": "1"} )
  if result==1:
    print( "Run terminated successfully after " + str(test.runtime) + "s ... ok" )
  else:
    print( "Run failed: " + test.last_output )
  serial_runtime = test.runtime
  
  for p in range(2,28,2):
    print( "Test for " + str(p) + " cores" )
    result = test.run( arguments, {"OMP_NUM_THREADS": str(p)} )
    speedup = serial_runtime/test.runtime
    if result==1 and speedup>1:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. with speedup of " + str(speedup) + " ... ok  (but you might want to tune it further)" )
    elif result==1:
      print( "Run terminated successfully after " + str(test.runtime) + "s, i.e. with deterioriated speedup ... check (but check runtimes)" )
    else:
      print( "Run failed: " + test.last_output )


def step5(zip):
  f = zip.open("report.pdf")
  content = f.read()
  f = open('report.pdf', 'wb')
  f.write(content)
  f.close()
  print( "========" )
  print( " Step 5" )
  print( "========" )
  print( "Extracted report.pdf ... ok" )



if __name__=="__main__":
  parser = argparse.ArgumentParser(description="""Parallel Scientific Programming I - 

 This is a small Python script to validate that the format of a submission
 is valid. It does not check for any correctness (though it gives some clues
 on the performance). If you use the script on Hamilton - which is the type of
 machine
 I plan to use to assess the submissions - you have to load appropriate modules.
 I use python/3.6.8 and then need an Intel compiler. The environment has recently
 been updated, so use the latest Intel compiler to get best performance:

 module load python/3.6.8 intel/2020.4;

 As always, I recommend to use a compute node rather than the login nodes. To do 
 this, you have to call something like

 salloc -N 1 -p test.q python3 validate.py test.zip;


""")
  parser.add_argument("--validate",    dest="validation_run",    action="store_true",      help="Validate whether submission is valid", default=True )
  parser.add_argument("zipfile",       help="zip archive" )
  args = parser.parse_args()

  try:
    with ZipFile(args.zipfile, 'r') as zip:
      #
      # 
      #
      print( "Open zip file " + args.zipfile )
      submitted_number_of_zip_files = len(zip.infolist())
      correct_number_of_zip_files   = 5
      if submitted_number_of_zip_files==correct_number_of_zip_files:
        print( "Number of archived files: " + str(submitted_number_of_zip_files) + " ... ok" )
      else:
        print( "Number of archived files: " + str(submitted_number_of_zip_files) + " ... wrong. Should be exactly " + str(correct_number_of_zip_files ) )

      step1(zip)        
      step2(zip)        
      step3(zip)
      step4(zip)
      step5(zip)
      
      print( """
      
Disclaimer: This is a sole sanity check, i.e. I'll run further tests on correctness and scalability 
when I mark the coursework. But the script ensures that your submission is formally correct plus it
also does some very basic checks.
 
""")
              
      arguments = "-fopenmp -O3 -xhost --std=c++0x"
      
      #
      # Clean-up
      #
      output_files = [ f for f in os.listdir(".") if "extracted" in f or ".cpp" in f]
      for f in output_files:
        os.remove(f)
      print( "Cleaned up all extracted files ... ok" )
  except Exception as e:
    print( "Sorry, could not validate/read input file: " + str(e) ) 

