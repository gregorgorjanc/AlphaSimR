#!/usr/bin/env Rscript

# Find the tskitr library directory and library file
# and return appropriate linker flags for use in Makevars
findTskitrLib <- function() {
  libdir <- system.file("libs", package = "tskitr")
  if (!nzchar(libdir)) {
    # nzchar() checks for non-zero character length
    stop("Unable to locate the tskitr shared library directory!")
  }

  if (.Platform$OS.type == "unix") {
    candidates <- c(
      "tskitr.so", # Unix/Linux
      "tskitr.dylib" # macOS
    )
  } else if (.Platform$OS.type == "windows") {
    candidates <- c(
      "tskitr.dll.a", # MinGW/Rtools (default)
      "tskitr.lib", # MSVC (backup to MinGW/Rtools)
      "tskitr.dll" # DLL (backup to MinGW/Rtools)
    )
  } else {
    stop("Unknown .Platform$OS.type!")
  }
  libpaths <- file.path(libdir, candidates)
  libfile <- libpaths[file.exists(libpaths)][1]
  if (length(libfile) < 1) {
    stop("Unable to locate the tskitr library file in ", libdir, "!")
  }

  # TODO: there were lots of trials to make the below work on macOS, so I
  #       assume we will need Linux & Windows specific flags as well!?
  # Embed rpath so the runtime loader can find tskitr.so.
  ret <- sprintf("-Wl,-rpath,%s %s", shQuote(libdir), shQuote(libfile))
  return(ret)
}

# Render a Makevars file from a template by replacing placeholders.
# @param template character path to the template Makevars file
# @param output character path to the output Makevars file
renderMakevars <- function(template, output) {
  if (!file.exists(template)) {
    stop("Template file does not exist: ", template, "!")
  }
  tskitrLib <- findTskitrLib()
  lines <- readLines(con = template)
  lines <- gsub(
    x = lines,
    pattern = "@TSKITR_LIB@",
    replacement = tskitrLib,
    fixed = TRUE
  )
  writeLines(text = lines, con = output)
  invisible(TRUE)
}

# renderMakevars(template = "this_should_fail", output = "before_getting_to_output")
if (.Platform$OS.type == "unix") {
  # readLines(con = "src/Makevars.in")
  success <- renderMakevars(
    template = "src/Makevars.in",
    output = "src/Makevars"
  )
} else {
  # readLines(con = "src/Makevars.win.in")
  success <- renderMakevars(
    template = "src/Makevars.win.in",
    output = "src/Makevars.win"
  )
}
if (!success) {
  stop("renderMakevars() failed!")
}
