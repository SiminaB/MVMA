# Multiple plot function - from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##make ARMA(1,1) correlation matrix
ARMAcor <- function(phi, rho, n)
{
  C <- matrix(1, nrow=n, ncol=n)
  for(i in 1:n)
  {
    for(j in 1:n)
    {
      if(i != j)
      {
        C[i,j] <- phi*rho^abs(i-j)
      }
    }
  }
  C
}

##get 1-denominator of efficiency for p=2 (note it's r2 for the 1st component, r1 for the 2nd component)
oneMinusDenomEff2 <- function(rho112, rho212, r2)
{
  rho112^2*r2+rho212^2*(1-r2)
}

##get 1-denominator of efficiency for p=3 (note it's r2 and r3 for the 1st component, need to permute them for
##other 2 components)
oneMinusDenomEff3 <- function(rho112, rho212, rho113, rho213, rho123, rho223, r2, r3)
{
  oneMinusDenomEff2(rho112, rho212, r2)+oneMinusDenomEff2(rho113, rho213, r3)+
    rho123^2*r2*r3+rho223^2*(1-r2)*(1-r3)-
    rho113^2*rho212^2*r3*(1-r2)-rho112^2*rho213^2*r2*(1-r3)-
    2*rho112*rho113*rho123*r2*r3-2*rho212*rho213*rho223*(1-r2)*(1-r3)+
    2*(rho112*rho113-rho123)*(rho212*rho213-rho223)*sqrt(r2*r3*(1-r2)*(1-r3))
}

##function to get correlation between X1j+X2j and X1j+X2j for any p
corSum <- function(rho112, rho212, r1, r2)
{
  rho112*sqrt(r1*r2) + rho212*sqrt((1-r1)*(1-r2))
}

##calculate efficiency in terms of rho112, rho212, r1, and r2
##for p=2
effCalc2 <- function(rho112, rho212, r1, r2)
{
  n1 <- 1-oneMinusDenomEff2(rho112, rho212, r2)
  n2 <- 1-oneMinusDenomEff2(rho112, rho212, r1)
  d <- 1-corSum(rho112, rho212, r1, r2)^2
  c(n1/d, n2/d)
}

##calculate efficiency in terms of rho112, rho212, etc
##for p=3
effCalc3 <- function(rho112, rho212, rho113, rho213, rho123, rho223, r1, r2, r3)
{
  n1 <- 1-oneMinusDenomEff3(rho112, rho212, rho113, rho213, rho123, rho223, r2, r3)
  n2 <- 1-oneMinusDenomEff3(rho112, rho212, rho123, rho223, rho113, rho213, r1, r3)
  n3 <- 1-oneMinusDenomEff3(rho113, rho213, rho123, rho223, rho112, rho212, r1, r2)

  d <- 1-corSum(rho112, rho212, r1, r2)^2-corSum(rho113, rho213, r1, r3)^2-
    corSum(rho123, rho223, r2, r3)^2+
    2*corSum(rho112, rho212, r1, r2)*corSum(rho113, rho213, r1, r3)*corSum(rho123, rho223, r2, r3)
  c(n1/d, n2/d, n3/d)
}

##create 2x2 var-cov matrix
varMat2 <- function(varX, varY, corXY)
{
  V <- matrix(0,2,2)
  V[1,1] <- varX
  V[2,2] <- varY
  V[1,2] <- V[2,1] <- corXY * sqrt(varX*varY)
  V
}

##create 3x3 var-cov matrix
varMat3 <- function(varX, varY, varZ, corXY, corXZ, corYZ)
{
  V <- matrix(0,3,3)
  V[1,1] <- varX
  V[2,2] <- varY
  V[3,3] <- varZ
  V[1,2] <- V[2,1] <- corXY * sqrt(varX*varY)
  V[1,3] <- V[3,1] <- corXZ * sqrt(varX*varZ)
  V[2,3] <- V[3,2] <- corYZ * sqrt(varY*varZ)
  V
}

##invert 2 by 2 var-cov matrix
invert2By2Var <- function(v)
{
  v11 <- v[1,1]
  v22 <- v[2,2]
  v12 <- v[1,2]
  d <- v11*v22-v12^2

  dInv <- matrix(0, 2, 2)
  dInv[1,1] <- v22
  dInv[2,2] <- v11
  dInv[1,2] <- dInv[2,1] <- -v12

  dInv <- 1/d*dInv

  dInv
}

##make panels figures S4-S5
panelFigS4S5 <- function(I, subsetSigma2, RelEffSubs) 
{
  ggplot(RelEffSubs[as.character(RelEffSubs$I) == as.character(I) &
                      as.character(RelEffSubs$Sigma2) == subsetSigma2,], 
         aes(x=p, y=RelEff))+
    geom_point(size=3.0, aes(color=randVar, shape=Estimate))+
    theme_bw(base_size = 20)+
    theme(axis.line = element_line(colour = "black"),
          plot.title = element_text(size = 15, hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.text.align = 0,
          ##scale_x_continuous(breaks=c(5,10,15,20)) +
          legend.position = c(0.15,0.2),
          axis.line.x = element_line(color="black", size = 0.5), ##this is to show axes - bug in this version of ggplot2
          axis.line.y = element_line(color="black", size = 0.5)) +
    guides(color=FALSE) + ##remove legend corresponding to color aesthetic
    ##scale_color_discrete(name = "Random scenario") +
    scale_shape_discrete(name ="",
                         labels = c(expression(paste(RelEff)), 
                                    expression(paste(RelEff^T)))) +
    ylim(0.73, 1) 
}