#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define BUFFER_SIZE 5
#define NUM_ITEMS 20

int buffer[BUFFER_SIZE];
int count = 0; // Number of items in the buffer
int in = 0;    // Index to insert an item into the buffer
int out = 0;   // Index to remove an item from the buffer

// Function to simulate the production of an item
int produce_item() {
    return rand() % 100; // Generating a random item for simplicity
}

// Function to simulate the consumption of an item
void consume_item(int item) {
    printf("Consumed item: %d\n", item);
}

// Producer function
void producer() {
    for (int i = 0; i < NUM_ITEMS; i++) {
        int item = produce_item();

        // Wait if the buffer is full
       
        while (count == BUFFER_SIZE);

        // Produce item
        buffer[in] = item;
        in = (in + 1) % BUFFER_SIZE;

        
        count++;

        printf("Produced item: %d\n", item);
    }
}

// Consumer function
void consumer() {
    for (int i = 0; i < NUM_ITEMS; i++) {
        // Wait if the buffer is empty
        #pragma omp critical
        while (count == 0);

        // Consume item
        int item = buffer[out];
        out = (out + 1) % BUFFER_SIZE;

        #pragma omp critical
        count--;

        consume_item(item);
    }
}

int main() {
    // Seed for random number generation
    srand(123);

    omp_set_num_threads(2);

    #pragma omp parallel sections
    {
       //add required pragma
	#pragma omp section
        {
            // Producer runs in one section
            producer();
        }

         //add required pragma
	#pragma omp section
        {
            // Consumer runs in another section
            consumer();
        }
    }

    return 0;
}

