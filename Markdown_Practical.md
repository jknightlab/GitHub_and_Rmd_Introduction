# R Markdown Practical
Katie Burnham  
7 February 2017  

In the Git practical, we used markdown for the README file. You can already see
that it is quite useful to make simple text documents, though mainly if you are
using GitHub. However, there are other benefits. 

## Why not just use Word?

1. Open R Studio again and add the following to your README file:

```
![](https://ipetcompanion.com/feedapuppy/styles/media/puppy.jpg)
```
Have you ever tried to format images in Word, and make them all stay where you
want to? We are currently working with a pretty simple document, but once you
are over a couple of pages, markdown makes this much more straightforwad. The 
rule is that you have the ! and square brackets, then the path to your image goes
in the normal brackets. This can be a local file, e.g.

```
![](my_image.jpg)
```

2. In the drop down menu next to the Preview button, click Preview in Word, or 
pdf. Then look at the different outputs you can get. Web pages and PDFs look
nice, right? But html and LaTeX are too complicated to learn just for a short
document. Well, markdown is a simpler alternative. As well as looking nice on
GitHub, you can email the file generated and email it to someone (e.g. Julian).
You can imagine that this is a nice way of sharing a summary of what you have been working on. 



















In your R Console (the bottom left panel), enter the following commands to 
install a plotting package that we will now use:

```
install.packages("ggplot2")
```

The code below simulates genotyping and gene expression data, which we can then 
test for an eQTL


```r
## Simulate 10 SNPs with different minor allele frequencies for 300 individuals
set.seed(100)
n <- 300
maf <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
prob <- cbind((1-maf)^2, 2*maf*(1-maf), maf^2)
geno <- apply(prob, 1, function(p) sample(0:2, size=n, prob=p, replace=TRUE))
geno <- data.frame(paste("sample", 1:n, sep="_"), geno)
names(geno) <- c("sample", paste("snp", 1:(ncol(geno)-1), sep="_"))

## Simulate gene expression for one gene per SNP.
## Expression is base + snp + noise
set.seed(101)
snp <- geno
errorMeans <- rnorm(nrow(snp))
expr <- sapply(snp[-1], function(g) 7 + 1.5*g) + 
		matrix(rnorm(nrow(snp)*(ncol(snp)-1), mean=errorMeans, sd=0.5), nrow=nrow(snp))
expr <- data.frame(snp$sample, expr)
names(expr) <- c("sample", paste("gene", 1:(ncol(expr)-1), sep="_"))
```

## Plotting gene expression by genotyping

A convenient way to display gene expression values by genotype is as box plots. These provide a good, non-parametric, indication of the distributions. To convey a sense of the frequency of each genotype in the sample it is useful to also add points for each individual to the plot. Below is an example of how this might look for each of the ten SNP/gene pairs.


```r
library(ggplot2)
genoLong <- tidyr::gather(geno, snp, genotype, -sample)
exprLong <- tidyr::gather(expr, gene, expression, -sample)
dataLong <- cbind(genoLong, exprLong["expression"])
dataLong$genotype <- as.factor(dataLong$genotype) 
ggplot(dataLong, aes(genotype, expression)) +
        geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
        geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + 
        facet_wrap(~snp) + theme_bw()
```

![](Markdown_Practical_files/figure-html/boxplots-1.png)<!-- -->

## Estimating SNP effects

To obtain estimates of the genotypic contribution to gene expression we fit a simple linear regression model of the form $E_i = \beta_0 + \beta G_i + \varepsilon$, where $E_i$ is the vector of gene expression values for gene $i$ and $G_i$ is the genotype vector for SNP $i$. We are interested in the estimate for $\beta$ which indicates the change in gene expression for each copy of the second allele.


```r
fit <- mapply(function(e, g) lm(e ~ g), expr[-1], geno[-1], SIMPLIFY=FALSE)
betaHat <- sapply(fit, coef)[2,]
betaHat 
```

```
##    gene_1    gene_2    gene_3    gene_4    gene_5    gene_6    gene_7 
## 0.9416734 1.3518024 1.1324493 1.5005237 1.4329068 1.2928385 1.5106722 
##    gene_8    gene_9   gene_10 
## 1.6050215 1.6038507 1.6162730
```

We use the function confint to obtain 95% confidence intervals of the estimated SNP effects.


```r
ci <- sapply(fit, confint, "g")
rownames(ci) <- c("lower", "upper")
ci
```

```
##           gene_1    gene_2    gene_3   gene_4   gene_5   gene_6   gene_7
## lower -0.1043169 0.8457378 0.5723167 1.156449 1.017273 1.026215 1.286452
## upper  1.9876638 1.8578670 1.6925819 1.844598 1.848541 1.559462 1.734892
##         gene_8   gene_9  gene_10
## lower 1.412641 1.426662 1.443334
## upper 1.797402 1.781040 1.789212
```

## Plotting results


```r
estimates <- data.frame(estimate=betaHat, t(ci), maf=maf)
ggplot(estimates, aes(x=maf)) + geom_hline(yintercept=1.5) + 
        geom_hline(yintercept=0, linetype="longdash") + 
        geom_errorbar(aes(ymin=lower, ymax=upper)) +
        geom_point(aes(y=estimate))  + theme_bw()
```

![](Markdown_Practical_files/figure-html/plots-1.png)<!-- -->

In this example all resulting confidence intervals include the true value^[although sometimes only just] but intervals for small minor allele frequencies are large (and in one case this means that 0 is included in the CI). As one would expect the uncertainty in the estimate, as measured by the length of the confidence interval, decreases with increasing minor allele frequency. However, even at high MAF considerable uncertainty remains and point estimates are somewhat lacking in accuracy, overestimating the true effect.
