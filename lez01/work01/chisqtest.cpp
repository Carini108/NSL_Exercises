/****************************************************************************
01.1.3 plot chi^2 (plotted later in python) 
*****************************************************************************/
   int M_blocks=100;
   vector<double> chi_2(M_blocks);
   int N_generations=pow(10,4);

   //calculate the chi^2 for increasingly thin blocks (sub-intervals)
   for (int number_of_blocks=1; number_of_blocks<=M_blocks; number_of_blocks++) {
      vector<double> frequencies(number_of_blocks);//bin populations
      for (int i=0; i<N_generations; i++){ //draw N_generations random numbers
         double x = rnd.Rannyu(0.,1.);
         for (int j=0; j<number_of_blocks; j++) { //this cycle checks whether the random number falls within [j/n,(j+1)/n)
            if ( x>=(static_cast<double>(j)/static_cast<double>(number_of_blocks)) && x<(static_cast<double>((j+1))/static_cast<double>(number_of_blocks)) ) {
               frequencies[j]+=1;//increases the population of the corresponding bin
               break;//stops the for cycle, since x can fall in one bin only
            }
         }
      }
      chi_2[number_of_blocks - 1] = 0; // Corrected index! vectors are filled from position 0
      for (int k = 0; k < number_of_blocks; k++)
        chi_2[number_of_blocks - 1] += pow(frequencies[k] - static_cast<double>(N_generations) / static_cast<double>(number_of_blocks), 2) / (static_cast<double>(N_generations) / static_cast<double>(number_of_blocks));
      cout << chi_2[number_of_blocks - 1] << endl; // Corrected index!
   }
   Print( chi_2 , "chi_squares.dat" );
   chi_2.clear();

/*
## 01.1.3: 

Suppose $n$ numbers are <font color="darkorange">uniformly sampled</font> on an interval $[a,b]$. Then $np$ of those numbers are <font color="darkorange">expected</font> to fall <font color="darkorange">within a sub-interval</font> $[a',b']\subseteq[a,b]$, where the probability $p$ is obtained as the ratio between the measures of the intervals, i.e. $p=(b'-a')/(b-a)$.

If, for example, $[0,1]$ is partitioned into $M$ identical sub-intervals, then $p=1/M$ and the <font color="darkorange">$\chi^2$ Pearson's cumulative test</font> may warn whether the hypothesis that the numbers are <font color="darkorange">not drawn</font> from a uniform distribution can be safely <font color="darkorange">discarded</font>. 

The plot shows the <font color="darkorange">increase of $\chi^2_j$</font> as the number $j$ of sub-intervals rises from $1$ (no partition is done) up to $M=100$, being
$$\chi^2_j = \sum_{i=1}^j \frac{\left( n_i - n/j \right)^2}{n/j}$$

with $n_i$ being the tally of numbers fallen within the $i$-th sub-interval and $n$ the number of draws (here $10^4$) at fixed $j$. On average, $(n_i - n/j)^2 \simeq n/j$, therefore $\chi^2_j \simeq j$, i.e. <font color="darkorange">the number of sub-intervals</font>.

For the sake of <font color="darkorange">comparison</font> of the obtained histogram with the actual distribution, see the <a href="https://en.wikipedia.org/wiki/Chi-squared_distribution"> Wikipedia page</a>.
*/