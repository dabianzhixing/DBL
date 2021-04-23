# DBL
This is the code for the paper "DBL: Efficient Reachability Queries on Dynamic Graphs".
We use Boost library to serialize the index. So you need to install Boost before using. The version we used is Boost_1_73_0.
Please modify the Makefile to guarantee the correct path to Boost file. Then, just enter "make" to compile the code.

After compiling, there will be three file: construct, query, and insert.
construct: construct the index and output the index file.
query: read the index file and query file, and execute the query.
insert: read the index file and insert file, and update the index.

All the input file(inculde dataset, query file, and insert file) are preprocessed into edge list. That is, every edge is denoted as the start node and the end node.
For example,
1 3
5 4
3 2
6 3
......

To construct the index
./construct ./path_to_the_dataset ./path_to_the_index

To execute the query
./query ./path_to_the_index ./path_to_the_queryfile

To execute edge insertion
./insert ./path_to_the_index ./path_to_the_insertionfile Number_of_Insertions

The default size of DL label and BL label is 64. You could tune the label size to achieve better performance. 
