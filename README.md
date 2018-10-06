# Queing-Theory
A simulation of a custom queing policy called donation policy.
This project is a done as an excersice for the network modeling course.

We simulate the queing policies FCFS (first come first serve) and SJF (shortest job first) under different load values, 
then we propose a new queing policy that we call donation policy and we simulate it.

The donation policy is a combination between the SJF policy and the processor sharing policy, it works as the following:

We assume we have 2 types of jobs BIG and SMALL ones each time we get a big job we do as SJF do, we through it to the back of
all the small jobs but this time we will let each small job donate thisbig job an amount X of time, This means that the server will serve 
the big job for an amount X then serve the small job...
this happens repeatidly so as the small jobs are served the big job is getting smaller and smaller...
the amount of donation X is an hyperparameter that i tuned to find an optimal value, where it appears to be 2.

Finally we evaluate each of the policies by drawing:
* A graph showing the time delay for a new job to get finished with respect to the load on the system.
* A histogram showing the frequancy of the number of jobs in the system after each iteration of each of the policies.
