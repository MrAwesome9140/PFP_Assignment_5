#include <pthread.h> /*used in other parts of the assignment */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>  /* for uint64  */
#include <time.h>    /* for clock_gettime */
#include <math.h>
#include <atomic>    /*used in other parts of the assignment */

#define NUM_POINTS 1000000000


// Integration Approximation
int NUM_THREADS = 1;

pthread_t* threads;
int* thread_ids;
double* sum;

double pi = 0.0;
pthread_mutex_t mutex;
pthread_barrier_t barrier;

// Prefix Sums
int NUMS_SIZE = 100000;

double* nums;
double* prefix_sums;

// Atomic Intructions
std::atomic<double> pi_atomic{0.0};

/**
 * @brief Adds the given value to the atomic variable pi_atomic.
 * 
 * @param bar The value to be added to pi_atomic.
 */
void add_to_pi(double bar) {
  auto current = pi_atomic.load();
  while (!pi_atomic.compare_exchange_weak(current, current + bar));
}

double p(double x) {
  return (6.0/sqrt(1-x*x));
}

double g(double x) {
  return sqrt(1-x*x);
}

// Computes the approximate value of pi/2 using sequential integration method.
// The function computes the approximate integral of sqrt(1-x^2) over [-1, 1] using NUM_POINTS.
// The execution time of the function is printed to the console.
// Returns 0.
double seqIntApprox() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  double step = 2.0 / NUM_POINTS;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);
 
  // Compute approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = -1;
  for (int i = 0; i < NUM_POINTS; i++) {
    x = step * i - 1;
    pi += step * sqrt(1 - (x * x));
  }
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
 
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
 
  printf("Sequential pi approximation elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int) execTime);
  printf("Variation: Sequential\n");
  printf("pi / 2 = %.20f\n\n", pi);

  return 0;
}

/**
 * This function updates the value of pi using the Monte Carlo method without synchronization.
 * It computes the approximate integral of sqrt(1-x^2) over [-1, 1] for a portion of the total number of points.
 * @param arg A pointer to an integer representing the thread ID.
 * @return NULL
 */
void* updatePi_No_Synch(void* arg) {
  double step = 2.0 / NUM_POINTS;

  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    x = i * step - 1;
    pi += step * sqrt(1-x*x); 
  }

  return NULL;
}

/**
 * This function updates the value of pi using mutex locks to ensure thread safety.
 * It computes the approximate integral of sqrt(1-x^2) over [-1, 1] for a specific thread.
 * @param arg A pointer to an integer representing the thread ID.
 * @return NULL
 */
void* updatePi_Mutex(void* arg) {
  double step = 2.0 / NUM_POINTS;
  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    x = i * step - 1;
    pthread_mutex_lock(&mutex);
    pi += step * sqrt(1-x*x); 
    pthread_mutex_unlock(&mutex);
  }

  return NULL;
}

/**
 * This function updates the value of pi using an atomic operation.
 * It computes the approximate integral of sqrt(1-x^2) over [-1, 1] for a portion of the total number of points.
 * @param arg A pointer to an integer representing the thread ID.
 * @return NULL
 */
void* updatePi_Atomic(void* arg) {
  double step = 2.0 / NUM_POINTS;
  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    x = i * step - 1;
    add_to_pi(step * sqrt(1-x*x));
  }

  return NULL;
}

/**
 * This function updates the global array "sum" with each thread's part of the approximate integral of sqrt(1-x^2) over [-1, 1].
 * @param arg A pointer to the thread ID.
 * @return NULL
 */
void* updatePi_GlobalArray(void* arg) {
  double step = 2.0 / NUM_POINTS;
  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    x = i * step - 1;
    sum[thread_id] += step * sqrt(1-x*x);
  }

  return NULL;
}

/**
 * This function updates the value of pi using a local sum.
 * It computes the approximate integral of sqrt(1-x^2) over [-1, 1] using multiple threads.
 * Each thread computes its own part of the integral and updates the corresponding element in the global sum array.
 * @param arg A pointer to an integer representing the thread ID.
 * @return NULL
 */
void* updatePi_LocalSum(void* arg) {
  double step = 2.0 / NUM_POINTS;
  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double x = 0.0;
  double local_sum = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    x = i * step - 1;
    local_sum += step * sqrt(1-x*x);
  }

  sum[thread_id] = local_sum;

  return NULL;
}

/**
 * This function updates the value of pi using local sum and barrier synchronization.
 * @param arg A pointer to the thread ID.
 * @return NULL
 */
void* updatePi_LocalSum_Barrier(void* arg) {
  double step = 2.0 / NUM_POINTS;

  int thread_id = *((int*) arg);

  // Compute this threads part of the approximate integral of sqrt(1-x^2) over [-1, 1]
  double local_sum = 0.0;
  for (int i = thread_id; i < NUM_POINTS; i += NUM_THREADS) {
    double x = i * step - 1;
    local_sum += step * sqrt(1 - (x * x));
  }

  sum[thread_id] = local_sum;
  pthread_barrier_wait(&barrier);

  return NULL;
}

/**
 * Calculates the approximation of pi in parallel using different synchronization methods.
 * 
 * @param synch_type an integer representing the synchronization method to be used:
 *                   1 - No Synchronization
 *                   2 - Mutex
 *                   3 - Atomic Instructions
 *                   4 - Global Array
 *                   5 - Local Sum
 *                   6 - Local Sum with Barrier
 */
void parallelIntApprox(int synch_type) {

  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  double step = 2.0 / NUM_POINTS;

  // Intialize synchronization variables
  pthread_mutex_init(&mutex, NULL);
  pthread_barrier_init(&barrier, NULL, NUM_THREADS + 1);

  void* (*updatePi)(void*);

  switch (synch_type) {
    case 1:
      updatePi = updatePi_No_Synch;
      break;
    case 2:
      updatePi = updatePi_Mutex;
      break;
    case 3:
      updatePi = updatePi_Atomic;
      break;
    case 4:
      updatePi = updatePi_GlobalArray;
      break;
    case 5:
      updatePi = updatePi_LocalSum;
      break;
    case 6:
      updatePi = updatePi_LocalSum_Barrier;
      break;
  }

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

  // Create threads
  for (int i = 0; i < NUM_THREADS; i++) {
    thread_ids[i] = i;
    pthread_create(&threads[i], NULL, updatePi, thread_ids + i);
  }

  if (synch_type == 6) {
    // Wait for all threads to finish
    pthread_barrier_wait(&barrier);
  } else {
    for (int i = 0; i < NUM_THREADS; i++) {
      pthread_join(threads[i], NULL);
    }
  }

  if (synch_type == 5 || synch_type == 6) {
    // Combine results from threads
    for (int i = 0; i < NUM_THREADS; i++) {
      pi += sum[i];
    }
  }
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;

  if (synch_type == 3) {
    pi = pi_atomic.load();
  }
 
  switch (synch_type) {
    case 1:
      printf("Variation: No Synchronization\n");
      break;
    case 2:
      printf("Variation: Mutex\n");
      break;
    case 3:
      printf("Variation: Atomic Instructions\n");
      break;
    case 4:
      printf("Variation: Global Array\n");
      break;
    case 5:
      printf("Variation: Local Sum\n");
      break;
    case 6:
      printf("Variation: Local Sum with Barrier\n");
      break;
  }
  printf("Number of threads = %d\n", NUM_THREADS);
  printf("Parallel pi approximation elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int) execTime);
  printf("pi / 2 = %.20f\n\n", pi);
}

/**
 * Computes prefix sums for a given section of an array using multiple threads.
 * @param arg A pointer to an integer representing the thread ID.
 * @return NULL.
 */
void* prefixSum(void* arg) {
  int thread_id = *((int*) arg);

  // Compute which segment this thread will compute prefix sums for
  int start_index = thread_id * (((double) NUMS_SIZE) / NUM_THREADS);
  int end_index = (thread_id + 1) * (((double) NUMS_SIZE) / NUM_THREADS);

  // Compute prefix sums for assigned section for this thread
  prefix_sums[start_index] = nums[start_index];
  for (int i = start_index + 1; i < end_index; i++) {
    prefix_sums[i] = prefix_sums[i - 1] + nums[i];
  }

  // Wait for all threads to complete their assigned section
  pthread_barrier_wait(&barrier);

  // Compute prefix sums for the sections 1 through NUM_THREADS - 1
  int segment_size = NUMS_SIZE / NUM_THREADS;
  for (int i = 1; i < NUM_THREADS; i++) {
    // Find this threads subsection of the segment
    start_index = (i * segment_size) + thread_id * (((double) segment_size) / NUM_THREADS);
    end_index = (i * segment_size) + (thread_id + 1) * (((double) segment_size) / NUM_THREADS);

    // Compute prefix sums for this subsection
    double left_sum = prefix_sums[(i * segment_size) - 1];
    for (int j = start_index; j < end_index; j++) {
      prefix_sums[j] += left_sum;
    }

    // Watit for all threads to complete their subsection
    pthread_barrier_wait(&barrier);
  }

  return NULL;
}

/**
 * Calculates the parallel prefix sums of an array using multiple threads.
 * Initializes synchronization variables, creates threads, waits for all threads to finish,
 * and calculates the elapsed process CPU time.
 */
void parallelPrefixSums() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  // Intialize synchronization variables
  pthread_barrier_init(&barrier, NULL, NUM_THREADS + 1);

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

  // Create threads
  for (int i = 0; i < NUM_THREADS; i++) {
    thread_ids[i] = i;
    pthread_create(&threads[i], NULL, prefixSum, &thread_ids[i]);
  }

  // Wait for all threads to finish
  for (int i = 0; i < NUM_THREADS; i++) {
    pthread_barrier_wait(&barrier);
  }
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
 
  printf("Number of threads = %d\n", NUM_THREADS);
  printf("Size of array = %d\n", NUMS_SIZE);
  printf("parallel prefix sums elapsed process CPU time = %llu nanoseconds\n\n", (long long unsigned int) execTime);
}

/**
 * Computes all prefix sums sequentially and prints the elapsed process CPU time in nanoseconds.
 * @return void
 */
void seqPrefixSums() {
  uint64_t execTime; /*time in nanoseconds */
  struct timespec tick, tock;

  clock_gettime(CLOCK_MONOTONIC_RAW, &tick);

  // Compute all prefix sums sequentially
  prefix_sums[0] = nums[0];
  for (int i = 1; i < NUMS_SIZE; i++) {
    prefix_sums[i] = prefix_sums[i - 1] + nums[i];
  }
 
  clock_gettime(CLOCK_MONOTONIC_RAW, &tock);
  execTime = 1000000000 * (tock.tv_sec - tick.tv_sec) + tock.tv_nsec - tick.tv_nsec;
 
  printf("sequential prefix sums elapsed process CPU time = %llu nanoseconds\n", (long long unsigned int) execTime);
}

int main(int argc, char *argv[]) {
  int synch_type = 6;

  // Integration Approximation
  if (argc == 2) {
    // Run all numbers of threads for this synchronization type
    synch_type = atoi(argv[1]);
    if (synch_type == 0) {
      // Run the sequential version
      seqIntApprox();
    } else {
      // Run the specified synchronization type with all numbers of threads
      for (NUM_THREADS = 1; NUM_THREADS <= 8; NUM_THREADS *= 2) {
        threads = (pthread_t*) malloc(NUM_THREADS * sizeof(pthread_t));
        thread_ids = (int*) malloc(NUM_THREADS * sizeof(int));
        sum = (double*) malloc(NUM_THREADS * sizeof(double));

        parallelIntApprox(synch_type);

        free(threads);
        free(thread_ids);
        free(sum);
      }
    }
  }
  else if (argc >= 3) {
    synch_type = atoi(argv[1]);
    NUM_THREADS = atoi(argv[2]);
  
    if (synch_type == -1) {
      // Run all synchronization types with all numbers of threads
      for (synch_type = 1; synch_type <= 6; synch_type++) {
        for (NUM_THREADS = 1; NUM_THREADS <= 8; NUM_THREADS *= 2) {
          threads = (pthread_t*) malloc(NUM_THREADS * sizeof(pthread_t));
          thread_ids = (int*) malloc(NUM_THREADS * sizeof(int));
          sum = (double*) malloc(NUM_THREADS * sizeof(double));
          parallelIntApprox(synch_type);
          free(threads);
          free(thread_ids);
          free(sum);
        }
      }
    } else if (synch_type == 0) {
      // Run the sequential version
      seqIntApprox();
    } else if (synch_type > 0 && synch_type <= 6) {
      // Run the specified synchronization type with the specified number of threads
      threads = (pthread_t*) malloc(NUM_THREADS * sizeof(pthread_t));
      thread_ids = (int*) malloc(NUM_THREADS * sizeof(int));
      sum = (double*) malloc(NUM_THREADS * sizeof(double));

      parallelIntApprox(synch_type);

      free(threads);
      free(thread_ids);
      free(sum);
    }
  }

  // Experiment to find optimal h
  // int h = 540830219;
  // while (abs(seqIntApprox(g, -1.0, 1.0, h) / M_PI) > 0.01)
  //   h += 1;
  // printf("h = %d\n", h);

  // Sequential Prefix Sums
  // int nums_sizes[4] = {100000, 500000, 1000000, 2000000};
  // for (int i = 0; i < 4; i++) {
  //   NUMS_SIZE = nums_sizes[i];
  //   nums = (double*) malloc(NUMS_SIZE * sizeof(double));
  //   prefix_sums = (double*) malloc(NUMS_SIZE * sizeof(double));

  //   for (int i = 0; i < NUMS_SIZE; i++) {
  //     nums[i] = (rand() / RAND_MAX) * NUMS_SIZE;
  //   }

  //   seqPrefixSums();
  //   free(nums);
  //   free(prefix_sums);
  // }


  // Prefix Sums
  if (argc == 4) {
    NUMS_SIZE = atoi(argv[3]);

    if (NUMS_SIZE > 0) {
      // Compute prefix sums with specified array size for all numbers of threads
      nums = (double*) malloc(NUMS_SIZE * sizeof(double));
      prefix_sums = (double*) malloc(NUMS_SIZE * sizeof(double));

      for (int i = 0; i < NUMS_SIZE; i++) {
        nums[i] = (rand() / RAND_MAX) * NUMS_SIZE;
      }

      for (NUM_THREADS = 1; NUM_THREADS <= 8; NUM_THREADS *= 2) {
        threads = (pthread_t*) malloc(NUM_THREADS * sizeof(pthread_t));
        thread_ids = (int*) malloc(NUM_THREADS * sizeof(int));
        parallelPrefixSums();
        free(threads);
        free(thread_ids);
      }

      free(nums);
      free(prefix_sums);
    } else if (NUMS_SIZE == 0) {
      // Compute prefix sums with all array sizes for all numbers of threads
      int nums_sizes[4] = {100000, 500000, 1000000, 2000000};
      for (int i = 0; i < 4; i++) {
        NUMS_SIZE = nums_sizes[i];
        nums = (double*) malloc(NUMS_SIZE * sizeof(double));
        prefix_sums = (double*) malloc(NUMS_SIZE * sizeof(double));

        for (NUM_THREADS = 1; NUM_THREADS <= 8; NUM_THREADS *= 2) {
          threads = (pthread_t*) malloc(NUM_THREADS * sizeof(pthread_t));
          thread_ids = (int*) malloc(NUM_THREADS * sizeof(int));

          for (int i = 0; i < NUMS_SIZE; i++) {
            nums[i] = (rand() / RAND_MAX) * NUMS_SIZE;
          }

          parallelPrefixSums();

          free(threads);
          free(thread_ids);
        }

        free(nums);
        free(prefix_sums);
      }
    }
  }

  return 0;
}