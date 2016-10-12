setwdNew <- function (dir) 
{
    if (!file.exists(dir)) {
        dir.create(dir, recursive = TRUE)
    }
    setwd(dir)
}




my.write.table <- function (values, file = file, head = "Identifier", row.names = TRUE, 
    col.names = TRUE, append = FALSE, quote = FALSE, sep = "\t", 
    na = "NA", digits = NA) 
{
    if (is.vector(values)) {
        values = as.data.frame(values)
    }
    if (!is.na(digits)) {
        if (is.data.frame(values)) {
            for (i in 1:length(values)) {
                if (is.numeric(values[[i]])) {
                  hasDecimals = suppressWarnings({
                    as.integer(values[[i]]) != values[[i]]
                  })
                  hasDecimals[is.na(hasDecimals)] = FALSE
                  if (any(hasDecimals)) {
                    values[[i]] = signif(values[[i]], digits = digits)
                  }
                }
            }
        }
        else {
            if (is.numeric(values)) {
                values = signif(values, digits = digits)
            }
        }
    }
    if (row.names) {
        if (col.names) {
            write.table(matrix(c(head, colnames(values)), nrow = 1), 
                file = file, sep = sep, quote = quote, col.names = FALSE, 
                row.names = FALSE, append = append, na = na)
            write.table(values, file = file, sep = sep, quote = quote, 
                col.names = FALSE, row.names = TRUE, append = TRUE, 
                na = na)
        }
        else {
            write.table(values, file = file, sep = sep, quote = quote, 
                col.names = FALSE, row.names = TRUE, append = append, 
                na = na)
        }
    }
    else {
        write.table(values, file = file, sep = sep, quote = quote, 
            col.names = col.names, row.names = FALSE, append = append, 
            na = na)
    }
}
