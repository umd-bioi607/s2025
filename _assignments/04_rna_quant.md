---
type: assignment
date: 2024-04-10T11:59:00+5:00
title: "Assignment #3: RNA-seq quantification"
published: false
#pdf: /static_files/assignments/asg.pdf
#attachment: /static_files/assignments/asg.zip
#solutions: /static_files/assignments/asg_solutions.pdf
#published: true
due_event: 
    type: due
    date: 2024-04-28T11:59:00+5:00
    description: 'Assignment #3 due'
---

You will implement the basic Expectation Maximization (EM) algorithm for transcript-level quantification from RNA-seq data.  You will be given the alignments in a custom format that is described in detail below.  The goal will be to implement the full EM algorithm (where each alignment has an associated latent variable).

I also **highly** recommend you refer to the [original RSEM paper](https://academic.oup.com/bioinformatics/article/26/4/493/243395) for more details about my the derivations below.  Also, [this post describing the RSEM model](https://ro-che.info/articles/2017-01-29-rsem) is quite good.

To get started, you can clone [this](https://github.com/umd-bioi607/assignment_3_sample_data) repository, which contains a `data` sub-directory with the data that is described below that you will use for this project.

## Assignment ##

The goal of this project is to implement the EM algorithm for transcript quantification from RNA-alignments.  

## Overall structure

You will submit your assignment as a tarball named `BIOI607_A3.tar.gz`.  When this tarball is expanded, it should create a
**single** folder named `BIOI607_A3`.  This folder must be created in the directory where the decompression (i.e. `tar xzvf`) is done, and must not be nested inside any other folders. The details of how you structure your "source tree" are up to you, but the following **must** hold (to enable proper automated testing of your programs).

 * There should be a file at the top level called `squant`.  This file should be executable, and I'd recommend making them executable Python scripts as in assignments 0 and 1. 
 
 * There should be a README.md file in the top level directory.  This README file should contain the following information.
     
     - What did you find to be the hardest part of this assignment?
     - What resources did you consult in working on this assignment (view this as a form of citation; you shouldn't _copy_ code directly from anywhere in your assignment, but if you consulted other sources please list them here).

**Turnin** : The assignment turnin will be handled using Gradescope.  

### The problem of transcript quantification ###

The problem of transcript quantification from RNA-seq data is a basic computational challenge in bioinformatics, as we discussed in class.  Among other issues, a main difficulty arises as sequencing reads may align equally well to multiple transcripts / locations.  Determining the abundance of different transcripts depends upon knowing how many reads arise from these transcripts during sequencing, while determining the probability that a given multi-aligning read should be assigned to a specific transcript depends upon knowing the abundance of the different transcripts in the sample.  This leads to a "chicken-and-egg" problem that is commonly tackled by the EM algorithm (or other statistical inference approaches).

### Your data ###

You are provided with a sample alignment file from a specific RNA-seq sample.  In order to simplify the non-essential parts of the project, you won't be required to parse the underlying SAM / BAM file.  Instead the alignments are provided in a simplified text-based format.  This file can be found in the `data` folder of the repository and it is called `alignments_sample.txt.gz` (**note**: this file is gzipped simply to make the size of the repository smaller.  You can unzip this file before giving it to your program if you don't wish to deal with a gzipped file directly).  The format of the file is as described below.  In the description, each element of the file is given a descriptive name followed by a `:` and its type.  For example `foo:int` means that `foo` is an `int`, not that the file contains a literal `:int` following `foo`.  Further, for the purposes of clarity ` <tab> ` means a tab character, not a tab surrounded by spaces.

```
number_of_transcripts:int
<transcript_record_1>
<transcript_record_2>
...
<transcript_record_m>
<alignment_block_for_read_1>
<alignment_block_for_read_2>
...
<alignment_block_for_read_n>
```

where a transcript record is of the form

```
transcript_name:string <tab> transcript_length:int
```

and an alignment block is of the form

```
number_of_alignments_for_read:int
aln1_txp:string <tab> aln1_ori:string <tab> aln1_pos:int <tab> aln1_alignment_prob:double
aln2_txp:string <tab> aln2_ori:string <tab> aln2_pos:int <tab> aln2_alignment_prob:double
...
alnk_txp:string <tab> alnk_ori:string <tab> alnk_pos:int <tab> alnk_alignment_prob:double
```

That is, the file has a header that lists each transcript's name and length, followed by an easy-to-parse block that encodes the alignments for each read.  For a given read, it can align to any number of places.
The first line of the alignment block gives the number of alignments for this read, and then for each alignment we simply provide:
    
* The name of the transcript on which the alignment occurs
* The orientation of the alignment (`f`|`rc`)
* The position of the leftmost base in the alignment 
* The conditional alignment probability (don't worry about how this is calculated for the time being)

### Your goal ###

**Your task is to write a basic read quantification tool, to determine the abundances of the transcripts listed in the header of this file based on the alignments themselves**.

### Your tool ###

You should implement a tool named `squant` (for simple quantifier).  The `squant` tool takes *2* mandatory arguments.  The first option is the path to a `.txt` format alignment file. The second is the path to the file where the output should be written.  For example:

~~~bash
$ squant alignments.txt quants.tsv
~~~

Your program should output the results in the following format.

```
name <tab> length <tab> est_frag
name_1:string <tab> length_1:double <tab> est_frag_1:double
...
name_m:string <tab> length_m:double <tab> est_frag_m:double
```

That is, the output is a simple 3-column TSV containing the transcript name, the length and the estimated number of reads arising from this transcript.

### The models ###

In your EM-algorithm, you are expected to model the transcript length and, the fragment length probability and the alignment probability.

Recall our main likelihood equation:

\\[ \mathcal{L}(\eta | \mathcal{F}) = \prod_{i=1}^n \sum_{t_j \in \text{alignments}(f_i)} \Pr(t_j \mid \eta ) \cdot \Pr(f_i \mid z_{ij} = 1) 
\\]

The first probability is simply given by our (estimated) transcript abundance i.e. \\(\Pr(t_j \mid \eta ) \propto\eta_j\\).  

The second probability \\(\Pr(f_i \mid z_{ij} = 1)\\) should, in your project, consist of 3 components.

 1. The probability of selecting a position (uniformly at random) along this fragment from which to draw the read
 2. The alignment probability (this is given directly in the file)
 3. The probability of drawing a fragment of the given length under the fragment length distribution (which, as you will see below, is a normal distribution with mean 200 and standard deviation 25).
 
In other words, \\(\Pr(f_i \mid z_{ij} = 1) = \Pr(p_j(f_i) \mid z_{ij} = 1) \cdot \Pr(a_j(f_i) \mid z_{ij} = 1) \cdot \Pr( \ell(f_i) \mid z_{ij} = 1) \\) where the first term here encodes the probability of selecting a position along the transcript to draw the fragment, second term here encodes the probability of the alignment itself (given as the 4th entry in each alignment record) and the third term is the probability of selecting a fragment compatible with the observed length of the one implied by this read.  The first probability \\(\Pr(p_j(f_i) \mid z_{ij} = 1)\\) is simply proportional to \\(\frac{1}{\ell_j}\\), since we assume a position is chosen uniformly at random.  The second probability here is trivial, since it is given.  Let's briefly cover the third. 

#### The fragment length distribution ###

First, we describe the fragment length distribution, since it will be essential to determine this probability.

Since this is a single-end experiment, we cannot derive the fragment length distribution directly from the data.  Therefore, I am telling you that it is generated from a normal distribution with mean \\(\mu = 200\\) and a standard deviation \\(\sigma = 25\\).  Recall that, even though we are observing some reads of _fixed length_ (100 in this case), they are just one end of a fragment that can be much longer.  The lengths of those underlying fragments is drawn according to this distribution.  We will actually assume that this distribution in _truncated_, in that there are no fragments of length < 0 or > 1000. We will use the notation \\(d\\) to represent the probability density function of our fragment length distribution.

#### The fragment length probability ###

We observe, for every alignment, one end of a paired end read.  In reality, the underlying fragment from which this read is drawn from one end, could have many different possible lengths.  Thus, we will consider them _all_, each appropriately weighted by its probability.

Assume that under a given alignment \\(a_j(f_i)\\) of fragment i to transcript j, the fragment aligns to the transcript in the forward orientation at position 500.  Further, assume that transcript \\(t_j\\) has a length of 1000 nucleotides.  What possible lengths could fragment i have had?  Well, it could have had any length from 0 to 500. Why, is the maximum possible length 500?  Because if the fragment were of length > 500, then the other end would be positioned beyond nucleotide 1000 of this transcript, but this transcript is only 1000 nucleotides long.  This means we could have had a fragment of length 0, **or** 1, **or**, 2, **or**, ..., **or** 500.  Recalling the rules for the probability of events A **or** B **or** C, etc. we can infer that the probability we are interested in is:

$$
\Pr(\ell(f_i) = 0) \text{ or } \Pr(\ell(f_i) = 1) \text{ or } \cdots \text{ or } \Pr(\ell(f_i) = 500) \\
= \Pr(\ell(f_i) = 0) + \Pr(\ell(f_i) = 1) + \cdots \text{ or } \Pr(\ell(f_i) = 500) \\
= \sum_{k=0}^{500} \Pr(\ell(f_i) = k) = \sum_{k=0}^{500} d(k)
$$

Another way to write this, if we simply take \\(D\\) to be the cumulative distribution function (CDF) of \\(d\\) is as

$$
 \Pr(\ell(f_i) \le 500) = \sum_{k=0}^{500} d(k) = D(k)
$$

This leads to a simple rule for determining \\(\Pr( \ell(f_i) \mid z_{ij} = 1, d)\\).

For a given alignment \\(a_j(f_i)\\), let \\(p_j(f_i)\\) be the position where the read alignment for read i begins on transcript j.  Then:

$$
\Pr( \ell(f_i) \mid z_{ij} = 1) = 
  \begin{cases}
  D(\ell(t_j) - p_j(f_i)) & \text{if} \text{ ori}(f_j) \text{ is f} \\\\
  D(p_j(f_i) + 100 & \text{if} \text{ ori}(f_j) \text{ is rc}
  \end{cases}
$$

Explaining this formula a bit, the first condition occurs when the read aligns in the forward orientation.  This is because, if the alignment for this read (the observed read) is on the forward strand, the situation is as below:

```
                -----[observed read]-->
---[transcript sequence]--------------------------------------------->
                                          <-----[unobserved mate]--
```

On the other hand, if the observed read is on the reverse complement strand, then the situation is as (again, forgive the programmer art):

```
                -----[unobserved mate]-->
---[transcript sequence]--------------------------------------------->
                                          <-----[observed read]--
```

So, if the observed alignment is to the forward strand, the fragment size is between the leftmost position of this read, and the rightmost position of the unobserved mate (which can go as far to the _right_ as the end of the transcript).  Conversely, if the observed alignment is to the reverse complement strand, the fragment size is between the rightmost end of this read and the leftmost position of the unobserved mate (which can go as far to the _left_ as the beginning of the transcript).

<!-- #### <a name="efflen"> Effective length </a> ####

The fact that the fragments have a non-zero length implies that, though each transcript has some real length in terms of nucleotides, it is _not_ true that a read can start at any place along a transcript.  For example, on a transcript of length 1000, there are very few fragments that could start at position 900.  This is because only fragments whose underlying fragment length are 100 or less could start at this position.  On the other hand, _many_ possible fragments could start at position 200, because if a read aligns (in the forward orientation) starting at position 200 it could have a fragment length of 0, 1, 2, ..., 800.  Any fragment length such that the other end does not exceed the transcript length is possible.  

To account for this fact, we define the effective length of each transcript based on the expected number of locations along this transcript that a fragment could start. Let \\(d\\) be the fragment length distribution. For sufficiently long transcripts (transcripts of true length > 1000), we can take the effective length to be \\(\hat{\ell}_i = \ell_i - E[d]\\).  For shorter transcripts, we have to be a bit more careful.  Specifically, we have to account that fragments longer than the transcript cannot be drawn from this transcript; we have to condition on the transcript length.  Let \\(d_k\\) be the distribution \\(d\\) conditioned on having a maximum fragment length of \\(k\\) (this is just taking the part of the distribution from 0 to k and re-normalizing it so it is a proper distribution).  Then, for any transcript of length k < 1000, we can define its effective length as \\(\hat{\ell}_i = \ell_i - E[d_k]\\). 

**Note**: You are required to write a function yourself to compute the effective length of each transcript from the given (actual) lengths and the provided parameters \\(\mu, \sigma\\) of the fragment length distribution.  However, to help you with debugging, I am providing the correct effective lengths of these transcripts.  Your values should be close to (e.g. within 1 or 2 of) these values.  These effective lengths are given in the file `elens.tsv` in the `data` folder.

#### A note on the simplified equivalence class model ####

Note that some of the probabilities we discussed above don't make sense in the equivalence class-based variant of the algorithm.  If we say all reads aligning to the same set of transcripts are _equivalent_ from the perspective of quantification, than each such read cannot have a distinct alignment probability or fragment length probability (the alignment score and position portions of the alignment record essentially go unused).

**Given** that this approach will be approximating the true likelihood we wish to evaluate anyway, _I encourage you to explore what choices you might make so this likelihood is as accurate as possible_.  For example, since all fragments in an equivalence class will have to use the same alignment probability, perhaps it doesn't make sense to include very unlikely alignment probabilities?  That is, if there is a read \\(r_i\\) that maps to transcripts 1, 2, 3, and 4 with conditional alignment probabilities 1.0, 1.0, 1.0, \\(1e^{-8}\\), perhaps that last alignment should be ignored for the purposes of defining the equivalence class.  I encourage you to play with these types of decisions and figure out what types of decision rules for equivalence class inclusion give you the best results. -->

## Assessing your results ##

In the `data` folder, along with `alignments.txt.gz`, you will find a file called `true_counts.tsv`.  This file provides the true number of reads coming from each transcript.  The goal of providing you with this file is to allow you to assess the accuracy of your quantification tool, and to determine it's accuracy as you work on your implementation.

Once you have generated your result file under each mode, you should compare your estimated counts with these true counts.  I encourage you to consider / explore different metrics for accuracy (you can find a bunch of different metrics implemented in [this simple python file](https://gist.githubusercontent.com/bshishov/5dc237f59f019b26145648e2124ca1c9/raw/c9ca9690c8e866c63680364ce14d387ad5e79c12/forecasting_metrics.py), or in standard packages like scikit learn).  

At the end of the day, your gradescope submission will be evaluated on the basis of the [spearman correlation](https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient) between the truth and prediction, the [mean absolute error](https://en.wikipedia.org/wiki/Mean_absolute_error) between the truth and prediction, and the [pearson correlation](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between the truth and prediction.
