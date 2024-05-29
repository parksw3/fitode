---
title: "reviewer comments for fitode paper"
---

## Reviewer 1

The manuscript "Fitting epidemic models to data -- a tutorial in memory of Fred Brauer" was a pleasure to read, and will make an excellent contribution to the literature. It is well written, appropriately balances breadth and depth, introduces a new computational tool, and it and fills a sizable information gap that exists right now for mathematicians (with limited or no background in applied statistics) -- and many statisticians -- who are interested in modern approaches to fitting nonlinear dynamic models to time series data. I commend the authors for this contribution, and for their remebrance of Fred.

I did not find any major faults with the paper, but have suggested some edits below that I feel should improve the paper. The first four comments are most substantive.

1. I urge the authors to be more explicit about the "observation process model"/"observation model" as an important link in connecting ODEs to data. In practice, data are often not direct observations of state variable, e.g., they might be the number of reported cases which might be proportional to the true case count, as already mentioned by the authors. It is therefore common practice, when building your likelihood function, to create an "observation model" of your data. This is typically done using a joint PDF/PMF parameterized by ODE model outputs, as this let's us define a likelihood function to use in an MLE framework. The authors first mention an "observation process model" in section 5, well after the framework has been introduced, and the authors' package (fitode) takes an observation model argument as illustrated in the code below line 357 and as described above equation (27). Earlier, explicit mention of the observation model as a concept would go a long way to
put those existing references into appropriate context.

I trust the authors judgement in making the above modification, but just in case it's helpful, here are my suggestions for how to more explicitly introduce the "observation model" concept: First, introduce the idea briefly with one or two sentences early in section 3. Something as simple as the following might suffice: "Here we assume our data are direct observations of one of our state variables. However, this is frequently not the case. When a more nuanced relationship defines the link between model and data, we can specify an \textbf{\textit{observation model}} that describes how our data values relate to the ODE trajectories. We will revisit this concept in the next section." Then, in section 4, mention it again in the context of constructing the likelihood function: Minor adjustments to the text preceding eq. (14) could reframe that derivation using the concept of an observation model more explicitly. With those two modifications in place, then the text leading up to
the "observation" option in in fitode (see eq. (27)) could be modified to refer to the observation model as an explicit concept.

**Done**.

2.  Notation: Somewhat related to observation model comments above is the use of "x()" for ODE state variable values, and "x[]" for data values. I find this to be atypical and probably not helpful to most readers.  I strongly encourage the authors to consider notation that is easier to grasp.  My recommendation: Using the upper case "X()" for data crossed my mind, however I think a better choice might be "y()", which greatly improves the readability compared to discerning "()" from "[]". While a lower level undergrad might find this use of "x vs y" notation a little problematic, I think the target audience here will find it much easier notation to follow than the current notation. Moreover, y is a commonly used symbol for the response values in statistics, so this use of "y()" will not be terribly foreign to readers familiar with existing statistical theory notation.

**Let's discuss this further. DJDE really likes this notation. BMB is not convinced.**

3. I'm confused by the $\Deltax[t_\ell}]$ term at the end of eq. (14), which doesn't quite follow from the assumptions in the preceding text. It seems as if the authors were at one point assuming that the probability of an observation x[t] can be approximated with a discrete distribution, letting $P(x=x[t]) = f(x[y], \sigma) \Delta x$ and then this discretized Normal distribution could be used to construct the likelihood function. The rest of the text, however, does not reflect this, nor does it otherwise explain what the $\Delta x$ term represents. If this was the intention, then there needs to be (1) clarification in the preceding text that this is what is being assumed, and (2) some follow-up text after eq. (17) clarifying how this \Delta x term turns into a constant in the log-likelihood function and can ultimately be ignored in the optimization step that follows. However, as it stands, I'm not convinced that this is the most practical way to introduce the concept of
likehood and MLE. Alternatively, please consider dropping the $\Delta x$ from eq. (14) and changing the left side of equation (14) to indicate that it's a joint density function, not a probability, and then define the likelihood from there with a followup sentence or footnote about how likelihoods are defined the discrete case usng PMFs instead of PDFs.  This approach (defining the likelihood using the joint PDF or PMF) is a common way to introduce the concept of a likelihood, and in this context, that approach might be preferable.

**Done** (not in a footnote; I think we should minimize our use of footnotes).

I should add that, for all three comments above, there is value in using notational conventions and other formalisms familiar to statisticians, where possible. This would help mathematicians who read this paper to more easily see how it relates to existing statistics literature, and likewise would make the paper more accessible to statisticians who have limited exposure to dynamic models and these techniques.


4. The discussion could include a bit more guidance regarding identifiability issues. For example, line 533, the parenthetical discription for "unidentifiable" might be reconsidered, and replaced with something that elaborates a bit more on the problem of (a) structural unidentifiability and (b) practical identifiability, each with their own description touching upon the fact that many ODE models (used in this context) are overparameterized and yield non-unique parameter estimates, even under ideal data assumptions. The second issue can arise for even strucuturally identifiable models which may be practically unidentifiable due to, e.g., a lack of data from certain parts of state space. A simple verbal example to illustrate this last point (if you wanted to go into that much detail) is to ask the reader to consider fitting the logistic growth model $dx/dt=rx(1-x/K)$, with known initial value $x_0$, to time series data that ends while the trajectory is only in the exponential
growth phase. In that case, the parameter $r$ will be confidently estimated, but $K$ will not since all sufficiently large $K$ values will appear to give equally good fits, which reflects that the data contain no information about that steady state value, $K$. Some or all of this might be better placed in the paragraph starting at line 600, where you discuss convergence issues. Consider also adding a brief mention of identifiability issues somewhere in section 4, and please consider adding some references for dealing with the different types of identifiability issues.

5. Related to the above comment, on line 603, it might be useful to call this procedure the "multistart method" (as it is sometimes called) and to mention that it is not only useful for diagnosing convergence issues (where you would see dissimilar "best" parameter sets with differing likelihood values), but it's also useful for detecting identifiability issues (here, differing best parameter sets would show very similar likelihood values).

**Done**

6. Consider mentioning some other distribution options above eq. (17), around line 242-244, like the various Generalized Poisson distributions. Include references (some of these are available in other R packages that might be worth mentioning). I find that mathematicians often are unaware that these other options exist, and so they only consider Normal, Poisson, and Negative Binomial.

**Added something. Not sure how important this really is ... all applications of more exotic count data distributions are for non-dynamic/regression models ...**

7. Minor inconsistencies in referring to equations: These include model (2), Equation (2), equation (2), distance (7), approximation (6), dR/dt (2), etc. I don't recall this journal having strict guidelines for equation referencing, so at a minimum I would ask the authors to consider dropping the instances of capital-E "Equation" to lower case.

**Not done**

8. There are a few places where table captions and R code run off into the margins.

**The bodies of Tables 1 and 2 were a bit too wide, I tweaked them a bit, not sure how much (more) effort this is worth ... Figure 1 caption is too long and overruns the page number/bottom margin**

9. Generally speaking, there might be room to include a few more references throughout, especially given the tutorial nature of this manuscript. Readers looking to apply this approach would benefit from some extra guidance finding literature related to some of the concepts or methods that are only briefly mentioned here.

**Any thoughts?**


## Reviewer 2 


This manuscript introduces the software package fitode, an R-based tool developed to aid the fitting of ODE models to observed time series data, particularly for epidemiological applications. The software serves as a practical response to the historical curiosity about how dynamic models fit to empirical data, a question posed by the late mathematical biologist Fred Brauer. The manuscript not only shows the functionality of fitode through examples involving compartmental epidemic models but also provides a tutorial approach to guide readers, presumably with a background in mathematics but less familiarity with statistical methods and optimization techniques. The discussion delves into the technical details of model fitting, parameter estimation, and the challenges inherent in this domain.

### Major Comments:

1. The authors acknowledge the issue of non-uniqueness in parameter values achieving the minimum error in model fitting. This is a critical limitation as it affects the reliability of the model predictions and the interpretability of the results. To address this challenge, the authors might consider incorporating regularization techniques which can help in constraining the parameter space and reducing the likelihood of overfitting. It would also be helpful if the manuscript could discuss strategies to identify and handle the impact of parameter correlations, which often contribute to this non-uniqueness.

**Partly done** (added para. on identifiability, expanded comment about regularization slightly)

2. While fitode is presented as a user-friendly tool tailored for epidemiologists, a systematic comparison with existing methods available in platforms like MATLAB or Berkeley Madonna is missing. Each of these tools has its strengths and limitations which could significantly influence user choice depending on their specific needs. For instance, MATLAB offers a broad range of built-in functions for optimization and model fitting along with high computational power, but it may not be as accessible due to licensing costs. On the other hand, Berkeley Madonna is known for its ease of use and speed in running complex dynamic models but might lack some of the advanced statistical tools provided by R. The authors should elaborate on the comparative advantages of fitode, possibly in terms of its ease of integration with other statistical methods in R, its specific utility for epidemiological models, or any unique features that address the nuances of fitting disease transmission
models to data.

**Not done. I don't think a "systematic comparison" is worth it, but we could try to say something brief. Looked at Berkeley Madonna; it does have a 'curve fitting' option, but there's almost nothing in the docs about what it's actually doing ...** See @Marc+2020, also https://github.com/topics/epidemiology?l=matlab

