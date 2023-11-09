Steps to run program:
1. run "make build" from command line
2. run the command "./pi_approx -1 -1" to run executable with all synchronization
   types and all numbers of threads

To run all experiments for pi approximation and prefix sums:
"./pi_approx -1 -1 0"
Note: This will take a decent amount of time, as the number of points is large so 
      execution time is slow for low thread counts and synchronization with true
      and false-sharing.

How to run different variations of pi_approx:

The executable is run with arguments in the form "./pi_approx s n p"

s - type of synchronization to run pi approximation with
Possible values of s and their corresponding meaning:
    0 - Sequential (No multi-threading)
    1 - No Synchronization
    2 - Mutex
    3 - Atomic Instructions
    4 - Global Array
    5 - Local Sum
    6 - Local Sum with Barrier

n - number of threads to run program with
    - If this value is left unspecified, the specified synchronization type will
      run with all thread counts

Ex: "./pi_approx 2 4" - This will run the mutex (2) variation of pi_approx with 4 threads
Ex 2: "./pi_approx 4" - This will run the global array (4) variation of pi_approx with
                        1, 2, 4, and 8 threads

p - size of array with which to run parallel prefix sums

Special case: 0 - Will run 100,000, 500,000, 1,000,000, and 2,000,000 sized double arrays
                  with parallel prefix sums, printing out execution times of all sizes
                  and thread counts.

Ex: "./pi_approx -1 -1 1000" - Will run all variations of pi approximation and will
                               run prefix sums with an array of 1000 doubles with 
                               thread counts of 1, 2, 4, and 8
Ex: "./pi_approx 5 8 0" - Will run local sum variation of pi approximation with 8 threads.
                          Will also run all combinations of array sizes and thread counts
                          for prefix sums.