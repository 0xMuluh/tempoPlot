#' Plot the origin of features over time
#'
#' \code{plotOrigin()} stacks, for each sample of a subject linked to a
#' partner, such as a recipient to the donor of a faecal transplant or an
#' infant to its mother, the shares of its features that were found in the
#' partner only, in the subject before an event only, in both, or in
#' neither (new).
#'
#' @param x a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'   object.
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
#' @param assay.type \code{Character scalar}. Specifies the assay of
#'   abundances. (Default: \code{"counts"})
#'
#' @param detection \code{Numeric scalar}. A feature is present in a
#'   sample when its value is above this. (Default: \code{0})
#'
#' @param value \code{Character scalar}. Stacks the shares of the present
#'   features (\code{"occupancy"}) or of their abundance
#'   (\code{"abundance"}). (Default: \code{"occupancy"})
#'
#' @param align.to \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} with the time of an event for each linked subject,
#'   such as the transplant, in the units of \code{time.col}. The samples
#'   of a subject at or before its event are its own reference, and time is
#'   drawn relative to the event. If \code{NULL}, subjects have no
#'   reference, and features are from the partner or new. Not available
#'   for visits. (Default: \code{NULL})
#'
#' @param layout \code{Character scalar}. Draws stacked bars (\code{"bar"})
#'   or stacked areas (\code{"area"}). (Default: \code{"bar"})
#'
#' @param level \code{Character scalar}. Plots each subject
#'   (\code{"subject"}) or the mean shares of each \code{group.by} group
#'   at each time (\code{"group"}). (Default: \code{"subject"})
#'
#' @param group.by \code{Character scalar}. Specifies the column from
#'   \code{colData(x)} that groups subjects when \code{level = "group"}. If
#'   \code{NULL}, all subjects form one group. (Default: \code{NULL})
#'
#' @param max.gap \code{Numeric scalar}. Breaks areas where consecutive
#'   times of a panel are further apart than this, in \code{time.col}
#'   units. For visits, \code{NULL} breaks them at visits without a sample;
#'   otherwise \code{NULL} never breaks them. (Default: \code{NULL})
#'
#' @param ... additional arguments passed to
#'   \code{\link[ggplot2:geom_bar]{geom_col}} or
#'   \code{\link[ggplot2:geom_ribbon]{geom_ribbon}}.
#'
#' @details
#' A feature of a sample is from the partner when it is present in any
#' sample of the partner, and from the subject (\code{"own"}) when it is
#' present in any reference sample of the subject; \code{"both"} when
#' both hold, and \code{"new"} when neither does. All samples of the
#' partner count, whatever their time; subset \code{x} to restrict them.
#' Samples without present features are left out.
#'
#' Strain sharing is computed from strain profiles upstream; with a strain
#' assay, the plot shows the origin of strains, as in Aggarwala et al.
#' (2021, doi:10.1038/s41564-021-00966-0, Fig. 3c), Schmidt et al. (2022,
#' doi:10.1038/s41591-022-01913-0, Fig. 2a) and Ferretti et al. (2018,
#' doi:10.1016/j.chom.2018.06.005, Fig. 2).
#'
#' @return
#' A \code{ggplot} object.
#'
#' @seealso \code{\link{plotPartner}} for the dissimilarity to the
#'   partner.
#'
#' @name plotOrigin
#'
#' @examples
#' # Illustration only: hitchip1006 has no linked design. Subjects are
#' # paired at random, so the shares show sharing by chance.
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
#'
#' # Mean shares of features from the partner at each visit, by sex
#' plotOrigin(tse, time.col = "visit", group = "subject",
#'     partner.col = "partner", level = "group", group.by = "sex")
#'
#' # Each subject's first visit as its reference, as areas
#' tse$visit.number <- as.integer(tse$visit)
#' tse$first <- 1L
#' plotOrigin(tse, time.col = "visit.number", group = "subject",
#'     partner.col = "partner", align.to = "first", layout = "area",
#'     level = "group", group.by = "sex")
NULL

#' @rdname plotOrigin
#' @importFrom SummarizedExperiment colData
#' @importFrom stats aggregate
#' @importFrom ggplot2 labs
#' @importFrom methods setMethod
#' @export
setMethod("plotOrigin", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, group, partner.col, assay.type = "counts",
            detection = 0, value = c("occupancy", "abundance"),
            align.to = NULL, layout = c("bar", "area"),
            level = c("subject", "group"), group.by = NULL, max.gap = NULL,
            ...) {
        value <- match.arg(value)
        layout <- match.arg(layout)
        level <- match.arg(level)
        partner <- .get_partner(x, group, partner.col)
        linked <- !is.na(partner)
        .check_longitudinal(x[, linked], time.col, group)
        .check_column(x, group.by)
        y <- .align_time(x[, linked], time.col, group, align.to)
        mat <- .get_assay(x, assay.type)
        present <- .present(mat, detection)
        w <- (if (value == "occupancy") present else mat * present)[,
            linked, drop = FALSE]

        origins <- if (is.null(align.to)) c("partner", "new") else
            c("partner", "both", "own", "new")
        subject <- as.character(colData(x)[[group]])
        ref <- !is.null(align.to) &
            .time_numeric(colData(y)[[time.col]]) <= 0
        own <- present[, linked, drop = FALSE]
        # Weight of each origin in each sample of y.
        sums <- matrix(0, length(origins), ncol(y))
        for (i in split(seq_len(ncol(y)), subject[linked])) {
            p <- rowSums(present[, subject == partner[linked][i[1L]],
                drop = FALSE]) > 0
            r <- rowSums(own[, i[ref[i]], drop = FALSE]) > 0
            from <- cbind(partner = p & !r, both = p & r, own = !p & r,
                new = !p & !r)[, origins, drop = FALSE]
            sums[, i] <- crossprod(from, w[, i, drop = FALSE])
        }
        keep <- !ref & colSums(sums) > 0
        shares <- prop.table(sums[, keep, drop = FALSE], 2L)

        cd <- as.data.frame(colData(y)[, unique(c(time.col, group,
            group.by)), drop = FALSE], optional = TRUE)
        df <- cd[rep(which(keep), each = length(origins)), , drop = FALSE]
        df$.feature <- factor(origins, levels = origins)
        df$.value <- as.vector(shares)
        panel <- group
        if (level == "group") {
            df <- aggregate(df[".value"], df[c(time.col, group.by,
                ".feature")], mean)
            panel <- group.by
        }
        if (layout == "area") {
            df <- .stack_composition(df, time.col, panel, layout, max.gap)
        }
        ylab <- if (value == "occupancy") "Share of present features" else
            "Share of abundance"
        .plot_composition(df, time.col, panel, layout, ylab, grey = "new",
            ...) +
            labs(x = paste0(time.col, if (!is.null(align.to))
                paste(" since", align.to)), fill = "origin")
    }
)
