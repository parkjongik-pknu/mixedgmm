# gmm_mixed


* **HOW TO INSTALL AND LOAD:**

```
# You can install this package directly from GitHub using `devtools`.

install.packages("devtools")

library(devtools)

pak::pak("parkjongik-pknu/mixedgmm")

library(mixedgmm)
```
<br>

* **HOW TO USE:**


```
gmm_mixed(X_Q = NULL, X_W = NULL, X_D = NULL, K = 3)

X_Q : continuous variables

X_W : categorical variables

X_D : discrete variables

You need to separate the variables into continuous, categorical, and discrete types, and scale them before fitting the model.
```

