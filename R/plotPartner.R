#' Plot the dissimilarity of subjects to their partners over time
#'
#' \code{plotPartner()} draws, for each sample of a subject linked to a
#' partner, such as an infant to its mother or a recipient to the donor of
#' a faecal transplant, its dissimilarity to the samples of its own partner
#' and to the samples of the other partners. Own partners closer than the
#' others show that partners share their microbiota.
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
#'
#' @param dissimilarity \code{Character scalar}. Specifies the
#'   sample by sample dissimilarity in \code{metadata(x)}, as stored by
#'   \code{\link[mia:getDissimilarity]{mia::addDissimilarity()}}.
#'
#' @param time.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that contains sampling times: numeric, \code{Date},
#'   \code{POSIXct} or a factor of visits. Only the times of linked
#'   subjects are used.
#'
#' @param group \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that identifies subjects.
#'
#' @param partner.col \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that gives the partner of each subject, a value of
#'   the \code{group} column, on all its samples, and \code{NA} on the
#'   samples of partners and of subjects without one.
#'
#' @param layout \code{Character scalar}. Draws one line per subject and
#'   comparison (\code{"line"}), or boxes (\code{"box"}) or violins
#'   (\code{"violin"}) per time point, own and other partners side by
#'   side. (Default: \code{"line"})
#'
#' @param show.individuals \code{Logical scalar}. When \code{TRUE}, draws
#'   subject lines and points; when \code{FALSE}, only boxes, violins or
#'   the \code{smooth} curve. (Default: \code{TRUE})
#'
#' @param smooth \code{Character scalar}. Specifies a summary of each
#'   comparison on top: \code{"none"}, \code{"mean"} or \code{"median"}
#'   (over the samples at each distinct time), or \code{"loess"} (loess
#'   curve). (Default: \code{"none"})
#'
#' @param facet.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} to facet by. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks a subject's line where
#'   consecutive samples are further apart than this, in \code{time.col}
#'   units. For visits, \code{NULL} breaks lines at visits without a
#'   sample; otherwise \code{NULL} never breaks lines. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to the main layer, as in
#'   \code{\link{plotTrend}}.
#'
#' @details
#' A subject's own partner value is its mean dissimilarity to all samples
#' of its partner; the other partners value is its mean dissimilarity to
#' all samples of the partners of other subjects, which is the expected
#' value under random pairing. Samples of partners do not need times.
#' Without other partners, as when one donor gives to all recipients,
#' only the own partner is drawn.
#'
#' Yassour et al. (2018, doi:10.1016/j.chom.2018.06.007, Fig. 1D) compare
#' the dissimilarity of infants to their own and to unrelated mothers;
#' Smillie et al. (2018, doi:10.1016/j.chom.2018.01.003, Fig. 2) follow
#' the similarity of recipients to their donors after transplantation.
#'
#' @return
#' A \code{ggplot} object.
#'
#' @seealso \code{\link{plotOrigin}} for the features shared with the
#'   partner.
#'
#' @name plotPartner
#'
#' @examples
#' # Illustration only: hitchip1006 has no linked design. Subjects are
#' # paired at random, so own and other partners should not differ.
#' data("hitchip1006", package = "miaTime")
#' tse <- hitchip1006
#' # Some profiles are stored twice; keep the first copy
#' counts <- SummarizedExperiment::assay(tse, "counts")
#' tse <- tse[, !duplicated(t(counts))]
#' tse <- tse[, tse$subject %in% names(which(table(tse$subject) >= 3))]
#' tse$visit <- factor(ave(tse$time, tse$subject, FUN = rank))
#' subjects <- unique(tse$subject)
#' half <- seq_len(length(subjects) %/% 2)
#' partner <- setNames(subjects[half], subjects[-half][half])
#' tse$partner <- unname(partner[as.character(tse$subject)])
#' tse <- mia::addDissimilarity(tse, method = "bray", assay.type = "counts")
#' plotPartner(tse, dissimilarity = "bray", time.col = "visit",
#'     group = "subject", partner.col = "partner", layout = "box")
NULL

#' @rdname plotPartner
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 labs
#' @importFrom methods setMethod
#' @export
setMethod("plotPartner", signature = c(x = "SummarizedExperiment"),
    function(x, dissimilarity, time.col, group, partner.col,
            layout = c("line", "box", "violin"), show.individuals = TRUE,
            smooth = c("none", "mean", "median", "loess"), facet.by = NULL,
            max.gap = NULL, ...) {
        layout <- match.arg(layout)
        smooth <- match.arg(smooth)
        partner <- .get_partner(x, group, partner.col)
        linked <- !is.na(partner)
        .check_longitudinal(x[, linked], time.col, group)
        .check_column(x, facet.by)
        .check_flag(show.individuals)
        d <- .get_dissimilarity(x, dissimilarity)

        subject <- as.character(colData(x)[[group]])
        partners <- unique(partner[linked])
        mean_to <- function(i, own) {
            to <- if (own) partner[i] else
                setdiff(partners, c(partner[i], subject[i]))
            mean(d[i, subject %in% to])
        }
        cd <- as.data.frame(colData(x)[, unique(c(time.col, group,
            facet.by)), drop = FALSE], optional = TRUE)
        cd$.i <- seq_len(ncol(x))
        cd <- .break_gaps(.order_time(cd[linked, , drop = FALSE], time.col,
            group), time.col, group, max.gap)
        df <- rbind(cd, cd)
        comparison <- c("own partner", "other partners")
        df$.partner <- factor(rep(comparison, each = nrow(cd)),
            levels = comparison)
        df$.y <- c(vapply(cd$.i, mean_to, numeric(1), own = TRUE),
            vapply(cd$.i, mean_to, numeric(1), own = FALSE))
        df <- df[!is.na(df$.y), , drop = FALSE]
        df$.line <- interaction(df[[group]], df$.segment, df$.partner,
            drop = TRUE)
        .plot_trend(df, time.col, layout, show.individuals, smooth, NULL,
            NULL, ".partner", facet.by, NULL, FALSE, ...) +
            labs(x = time.col, y = paste("Mean", dissimilarity),
                colour = "partner", fill = "partner")
    }
)

# Candidate for utils.R (shared by plotPartner and plotOrigin).
# The partner.col value of each sample, as character: a subject of group,
# constant within subjects, NA for subjects without a partner.
#' @importFrom SummarizedExperiment colData
.get_partner <- function(x, group, partner.col) {
    .require_column(x, group)
    .require_column(x, partner.col)
    subject <- as.character(colData(x)[[group]])
    partner <- as.character(colData(x)[[partner.col]])
    if (anyDuplicated(unique(data.frame(subject, partner))$subject)) {
        stop("'partner.col' must be constant within each subject.",
            call. = FALSE)
    }
    if (all(is.na(partner)) || !all(partner %in% c(subject, NA))) {
        stop("'partner.col' must give subjects of the 'group' column, ",
            "for at least one subject.", call. = FALSE)
    }
    partner
}
