###############################################################################
##                             ReacNorm R package                            ##
##                  Helpers functions throughout the package                 ##
##       ----------------------------------------------------------------    ##
##                           Pierre de Villemereuil                          ##
##       ----------------------------------------------------------------    ##
##                                     2024                                  ##
###############################################################################

## --------------------------------------------------------------- LICENCE ----

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

##  ---------------------------------------------------- Backend functions ----

## Non-Bessel-corrected variance
var_nocorrect <- function(vec) {
    N <- length(vec)

    ((N - 1) / N) * var(vec)
}

##  --------------------------------------------------- Frontend functions ----

## Utilitary function (exposed to the user) to format a shape function automatically
# Args: - expr: expression from which the gradient needs to the computed (expression)
#       - pars: the parameters of the reaction norm
# Value: A function that can be used to compute V_A using QGrn_va()
rn_generate_shape <- function(expr, pars) {
    # Generating the new parameter "names" for the function to provide QGglmm
    replacements <- stri_c("pars[",1:length(pars),", ]")
    patterns     <- stri_c("(?<![:alnum:])", pars, "(?![:alnum:])")

    # Computing the derivative
    expr2 <-
        expr |>
        as.character() |>
        stri_replace_all_regex(pattern = patterns,
                               replacement = replacements,
                               vectorise_all = FALSE)

    # Constructing the body of the function
    body <- parse(text = expr2)

    # Generating the list of arguments without defaults
    args <- list()
    args[["x"]] <- alist(x=)$x
    args[["pars"]] <- alist(pars=)$pars

    # Return the function to provide QGglmm
    eval(call("function", as.pairlist(args), body[[1]]), parent.frame())
}

## Utilitary function (exposed to the user) to compute a gradient automatically
# Args: - expr: expression from which the gradient needs to the computed (expression)
#       - dpars: the parameters wrt which we need to compute the derivative (character)
#       - allpars: list of all the parameters (character)
# Value: A function that can be used to compute V_A using QGrn_va()
rn_generate_gradient <- function(expr, dpars, allpars) {
    # Generating the new parameter "names" for the function to provide QGglmm
    replacements <- stri_c("pars[",1:length(allpars),", ]")
    patterns     <- stri_c("(?<![:alnum:])", allpars, "(?![:alnum:])")

    # Computing the derivative
    deriv <-
        lapply(dpars, \(p) D(expr, p)) |>
        as.character() |>
        stri_replace_all_regex(pattern = patterns,
                               replacement = replacements,
                               vectorise_all = FALSE)

    # Constructing the body of the function
    body <- parse(text = stri_c(c("matrix(rbind(",
                                  stri_c(deriv, collapse = ","),
                                  "), nrow = ",
                                  length(dpars),
                                  ", ncol = ncol(pars))"),
                                collapse = ""))

    # Generating the list of arguments without defaults
    args <- list()
    args[["x"]] <- alist(x=)$x
    args[["pars"]] <- alist(pars=)$pars

    # Return the function to provide QGglmm
    eval(call("function", as.pairlist(args), body[[1]]), parent.frame())
}

## Utilitary function (exposed to the user) to compute the second-order automatically
# Args: - expr: expression from which the gradient needs to the computed (expression)
#       - dpars: the parameters wrt which we need to compute the derivative (character)
#       - allpars: list of all the parameters (character)
# Value: A function that can be used to compute V_A using QGrn_va()
rn_generate_2diff <- function(expr, dpars, allpars) {
    # Generating the new parameter "names" for the function to provide QGglmm
    replacements <- stri_c("pars[",1:length(allpars),", ]")
    patterns     <- stri_c("(?<![:alnum:])", allpars, "(?![:alnum:])")

    # Computing the first derivative
    deriv <-
        lapply(dpars, \(p) D(expr, p))

    # Computing the second derivative and formatting
    d_deriv <-
        mapply(D, deriv, dpars) |>
        as.character() |>
        stri_replace_all_regex(pattern = patterns,
                               replacement = replacements,
                               vectorise_all = FALSE)

    # Constructing the body of the function
    body <- parse(text = stri_c(c("matrix(c(",
                                  stri_c(d_deriv, collapse = ","),
                                  "), nrow = ,",
                                  length(dpars),
                                   ", ncol = ncol(pars), ",
                                  ", byrow = TRUE)"),
                                collapse = ""))

    # Generating the list of arguments without defaults
    args <- list()
    args[["x"]] <- alist(x=)$x
    args[["pars"]] <- alist(pars=)$pars

    # Return the function to provide QGglmm
    eval(call("function", as.pairlist(args), body[[1]]), parent.frame())
}
