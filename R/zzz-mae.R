# MultiAssayExperiment methods. A plot of one experiment takes
# experiment (name or index) after the arguments of the
# SummarizedExperiment method, so positional calls mean the same for both,
# and calls that method on .get_experiment(x, experiment). plotSwimmer()
# and plotProcrustes() have their own methods, as they read several
# experiments. The file name sorts it after the plot files, so that each
# Rd shows the SummarizedExperiment method and its arguments first.

# Name of the experiment given by name or index in the
# MultiAssayExperiment x. arg names the argument in the error.
#' @importFrom MultiAssayExperiment experiments
.experiment_name <- function(x, experiment, arg = "experiment") {
    exps <- names(experiments(x))
    if (is.numeric(experiment) && length(experiment) == 1L &&
            experiment %in% seq_along(exps)) {
        experiment <- exps[experiment]
    }
    if (!(.is_one_name(experiment) && experiment %in% exps)) {
        stop("'", arg, "' must be the name or index of one experiment of ",
            "'x'.", call. = FALSE)
    }
    experiment
}

# One experiment of the MultiAssayExperiment x, whose samples get the rows
# of colData(x) of their primary names. Columns of the experiment's own
# colData are kept unless colData(x) has a column of the same name.
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment colData colData<-
.get_experiment <- function(x, experiment) {
    experiment <- .experiment_name(x, experiment)
    se <- experiments(x)[[experiment]]
    if (!inherits(se, "SummarizedExperiment")) {
        stop("'experiment' must be a SummarizedExperiment.", call. = FALSE)
    }
    cd <- colData(x)[.primary_names(x, experiment), , drop = FALSE]
    rownames(cd) <- colnames(se)
    own <- colData(se)
    colData(se) <- cbind(own[, setdiff(names(own), names(cd)), drop = FALSE],
        cd)
    se
}

# Method for generic f that plots one experiment.
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom methods setMethod
.mae_method <- function(f) {
    function(x, ..., experiment = 1L) {
        f(.get_experiment(x, experiment), ...)
    }
}

#' @rdname plotAlluvial
#' @template mae-experiment
#' @export
setMethod("plotAlluvial", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotAlluvial))

#' @rdname plotCalendar
#' @template mae-experiment
#' @export
setMethod("plotCalendar", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotCalendar))

#' @rdname plotComposition
#' @template mae-experiment
#' @export
setMethod("plotComposition", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotComposition))

#' @rdname plotDOC
#' @template mae-experiment
#' @export
setMethod("plotDOC", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotDOC))

#' @rdname plotDirection
#' @template mae-experiment
#' @export
setMethod("plotDirection", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotDirection))

#' @rdname plotFeatureScatter
#' @template mae-experiment
#' @export
setMethod("plotFeatureScatter", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotFeatureScatter))

#' @rdname plotFingerprint
#' @template mae-experiment
#' @export
setMethod("plotFingerprint", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotFingerprint))

#' @rdname plotGaps
#' @template mae-experiment
#' @export
setMethod("plotGaps", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotGaps))

#' @rdname plotHorizon
#' @template mae-experiment
#' @export
setMethod("plotHorizon", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotHorizon))

#' @rdname plotIndividuality
#' @template mae-experiment
#' @export
setMethod("plotIndividuality", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotIndividuality))

#' @rdname plotLagHeatmap
#' @template mae-experiment
#' @param target.experiment \code{Character scalar} or \code{integer
#'   scalar}. For a \code{MultiAssayExperiment}, the experiment holding
#'   \code{target}, looked up there as for a \code{SummarizedExperiment},
#'   in its first assay for a feature. Only samples whose primary sample is
#'   in both experiments are used. If \code{NULL}, \code{target} is
#'   looked up in \code{experiment}. (Default: \code{NULL})
#' @export
setMethod("plotLagHeatmap", signature = c(x = "MultiAssayExperiment"),
    function(x, target, ..., experiment = 1L, target.experiment = NULL) {
        experiment <- .experiment_name(x, experiment)
        se <- .get_experiment(x, experiment)
        if (!is.null(target.experiment)) {
            # The target becomes a colData column of the shared samples.
            name <- .experiment_name(x, target.experiment,
                "target.experiment")
            y <- .resolve_y(experiments(x)[[name]], target, NULL, "target")
            p1 <- .primary_names(x, experiment)
            p2 <- .primary_names(x, name)
            if (anyDuplicated(p2[p2 %in% p1])) {
                stop("'target.experiment' must have one sample per ",
                    "primary sample.", call. = FALSE)
            }
            se <- se[, p1 %in% p2]
            se[[target]] <- unname(y[match(p1[p1 %in% p2], p2)])
        }
        plotLagHeatmap(se, target, ...)
    }
)

#' @rdname plotLagProfile
#' @template mae-experiment
#' @export
setMethod("plotLagProfile", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotLagProfile))

#' @rdname plotLandscape
#' @template mae-experiment
#' @export
setMethod("plotLandscape", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotLandscape))

#' @rdname plotNetwork
#' @template mae-experiment
#' @export
setMethod("plotNetwork", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotNetwork))

#' @rdname plotPath
#' @template mae-experiment
#' @export
setMethod("plotPath", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotPath))

#' @rdname plotRanks
#' @template mae-experiment
#' @export
setMethod("plotRanks", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotRanks))

#' @rdname plotRecurrence
#' @template mae-experiment
#' @export
setMethod("plotRecurrence", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotRecurrence))

#' @rdname plotSpectrum
#' @template mae-experiment
#' @export
setMethod("plotSpectrum", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotSpectrum))

#' @rdname plotStats
#' @template mae-experiment
#' @export
setMethod("plotStats", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotStats))

#' @rdname plotSubjectScatter
#' @template mae-experiment
#' @export
setMethod("plotSubjectScatter", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotSubjectScatter))

#' @rdname plotSurvival
#' @template mae-experiment
#' @export
setMethod("plotSurvival", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotSurvival))

#' @rdname plotTBI
#' @template mae-experiment
#' @export
setMethod("plotTBI", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTBI))

#' @rdname plotTiles
#' @template mae-experiment
#' @export
setMethod("plotTiles", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTiles))

#' @rdname plotTimeEffect
#' @template mae-experiment
#' @export
setMethod("plotTimeEffect", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTimeEffect))

#' @rdname plotTimeLoadings
#' @template mae-experiment
#' @export
setMethod("plotTimeLoadings", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTimeLoadings))

#' @rdname plotTransitionMatrix
#' @template mae-experiment
#' @export
setMethod("plotTransitionMatrix", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTransitionMatrix))

#' @rdname plotTrend
#' @template mae-experiment
#' @export
setMethod("plotTrend", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTrend))

#' @rdname plotTurnover
#' @template mae-experiment
#' @export
setMethod("plotTurnover", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotTurnover))

#' @rdname plotVolcano
#' @template mae-experiment
#' @export
setMethod("plotVolcano", signature = c(x = "MultiAssayExperiment"),
    .mae_method(plotVolcano))
