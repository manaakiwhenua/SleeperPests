###############################################################################
### Native-R illustrative GIF: unmanaged 2-stage point invasion
###
### Uses only base R plus INApestPointTransitionMatrix.R.
### No terra, ggplot2, magick, gifski, Python or external GIF encoder required.
###
### Illustration:
###   * 3000 x 3000 continuous-space landscape
###   * 18 timesteps
###   * two biological stages
###   * stage 1 can move when transitioning to stage 2
###   * stage 2 produces new stage-1 propagules
###   * new propagules are highlighted only in their birth timestep
###   * source -> destination lines are drawn for reproductive dispersal
###   * source -> destination lines are drawn for stage-transition movement
###   * 1 second per animation frame
###############################################################################

find_project_root <- function() {
  candidates <- unique(normalizePath(c(getwd(), dirname(getwd())), winslash = "/", mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "models", "INApestPointTransitionMatrix.R")))
      return(candidate)
  }
  stop("Run this script from the extracted package root or its demos/ directory.")
}

project_root <- find_project_root()
output_dir <- file.path(project_root, "outputs")
demo_dir <- file.path(project_root, "demos")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(demo_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(project_root, "models", "INApestPointTransitionMatrix.R"))

###############################################################################
### 1. Model parameterisation
###############################################################################

# Transition matrix convention used by INApestPointTransitionMatrix:
# rows = target stage; columns = source stage.
# Here:
#   stage 1 has some stasis and frequent progression to stage 2;
#   stage 2 is persistent and reproductive.
A_demo <- matrix(c(
  0.25, 1.30,
  0.60, 0.82
), nrow = 2, byrow = TRUE)

# Full-suitability base-R grid. This demonstrates that the model can use the
# new native grid architecture even though habitat does not constrain this run.
habitat_demo <- INApestSpatialGrid(
  values = matrix(1, nrow = 60, ncol = 60),
  xmin = 0, xmax = 3000,
  ymin = 0, ymax = 3000
)

initial_points <- data.frame(
  x = c(1450, 1500),
  y = c(1525, 1500),
  stage = c(1L, 2L)
)

fit <- INApestPointTransitionMatrix(
  ModelName = 'Illustrative_TM_Unmanaged_2stage_nativeR',
  Nperm = 1L,
  Ntimesteps = 18L,
  Nstages = 2L,
  Weights = c(1, 1),
  Transition = A_demo,
  InitialPoints = initial_points,

  # Reproductive dispersal. Stage 2 is the reproductive source in this setup.
  SDDkernel = list(
    INApestPointKernelExponential(150),
    INApestPointKernelExponential(260)
  ),
  LDDkernel = NULL,
  LDDrate = 0,

  PropaguleEstablishment = 0.92,
  EnvEstabProb = 1,

  # Explicit movement on stage 1 -> stage 2 transition.
  TransitionKernels = list(
    INApestPointKernelExponential(90)
  ),
  TransitionHabitatSearch = FALSE,

  HabitatSuitability = habitat_demo,
  HabitatSearchRadius = 0,
  LocalK = Inf,
  KRadius = 0,

  # Unmanaged invasion.
  DetectionProb = c(0, 0),
  DetectionSD = 0,
  ManageProb = 0,
  ManageSD = 0,
  MortalityProb = c(0, 0),
  MortalitySD = 0,
  FecundityReduction = 0,
  FecundityReductionSD = 0,
  SpreadReduction = 0,
  SpreadReductionSD = 0,

  InfoTransferProb = 0,
  InfoRetentionProb = 1,
  InfoPersistenceSteps = Inf,

  DoProgress = FALSE,
  Seed = 20260818
)

###############################################################################
### 2. Save model outputs used by the animation
###############################################################################

write.csv(initial_points,
          file.path(output_dir, 'illustrative_unmanaged_initial_points_nativeR.csv'),
          row.names = FALSE)
write.csv(fit$PointHistory,
          file.path(output_dir, 'illustrative_unmanaged_point_history_nativeR.csv'),
          row.names = FALSE)
write.csv(fit$EventLog,
          file.path(output_dir, 'illustrative_unmanaged_event_log_nativeR.csv'),
          row.names = FALSE)
write.csv(fit$Summary,
          file.path(output_dir, 'illustrative_unmanaged_summary_nativeR.csv'),
          row.names = FALSE)

###############################################################################
### 3. Helpers for finding source locations from model output
###############################################################################

# Return the known location of an individual immediately before timestep t.
point_location_before <- function(id, timestep) {
  if (timestep <= 1L) {
    if (id >= 1L && id <= nrow(initial_points)) {
      return(c(initial_points$x[id], initial_points$y[id]))
    }
    return(c(NA_real_, NA_real_))
  }

  z <- fit$PointHistory[
    fit$PointHistory$id == id &
      fit$PointHistory$timestep == timestep - 1L,
    c('x', 'y'),
    drop = FALSE
  ]

  if (nrow(z) > 0L)
    return(c(z$x[nrow(z)], z$y[nrow(z)]))

  # Fallback for a point born during the immediately preceding sequence.
  ev <- fit$EventLog[
    fit$EventLog$point_id == id &
      fit$EventLog$timestep < timestep &
      fit$EventLog$event %in% c('establishment', 'external_incursion'),
    c('x', 'y', 'timestep'),
    drop = FALSE
  ]
  if (nrow(ev) > 0L) {
    ev <- ev[order(ev$timestep), , drop = FALSE]
    return(c(ev$x[nrow(ev)], ev$y[nrow(ev)]))
  }

  c(NA_real_, NA_real_)
}

# Return the source location for a reproductive establishment event.
# Prefer the parent's post-transition/current-timestep location when it exists;
# otherwise fall back to the previous-timestep location.
parent_location_for_reproduction <- function(parent_id, timestep) {
  if (is.na(parent_id))
    return(c(NA_real_, NA_real_))

  z <- fit$PointHistory[
    fit$PointHistory$id == parent_id &
      fit$PointHistory$timestep == timestep,
    c('x', 'y'),
    drop = FALSE
  ]
  if (nrow(z) > 0L)
    return(c(z$x[nrow(z)], z$y[nrow(z)]))

  point_location_before(parent_id, timestep)
}

###############################################################################
### 4. Tiny base-R raster drawing engine
###############################################################################

canvas_width <- 700L
canvas_height <- 700L
plot_margin <- 28L
xmin <- 0
xmax <- 3000
ymin <- 0
ymax <- 3000

# Colour indices used in the GIF palette.
COL_BACKGROUND <- 0L
COL_GRID <- 1L
COL_BORDER <- 2L
COL_STAGE1 <- 3L
COL_STAGE2 <- 4L
COL_PROPAGULE <- 5L
COL_REPRO_LINE <- 6L
COL_TRANSITION_LINE <- 7L
COL_FOUNDER <- 8L

# 256-entry global palette. Only the first entries are used.
palette_rgb <- matrix(255L, nrow = 256L, ncol = 3L)
palette_rgb[COL_BACKGROUND + 1L, ] <- c(255, 255, 255) # white
palette_rgb[COL_GRID + 1L, ] <- c(232, 232, 232)       # light grey
palette_rgb[COL_BORDER + 1L, ] <- c(75, 75, 75)        # dark grey
palette_rgb[COL_STAGE1 + 1L, ] <- c(72, 170, 205)      # blue/cyan
palette_rgb[COL_STAGE2 + 1L, ] <- c(0, 150, 136)       # teal
palette_rgb[COL_PROPAGULE + 1L, ] <- c(230, 75, 53)    # orange-red
palette_rgb[COL_REPRO_LINE + 1L, ] <- c(246, 166, 154) # pale coral
palette_rgb[COL_TRANSITION_LINE + 1L, ] <- c(183, 168, 150) # taupe
palette_rgb[COL_FOUNDER + 1L, ] <- c(20, 20, 20)       # near-black

xy_to_pixel <- function(x, y) {
  px <- plot_margin +
    round((x - xmin) / (xmax - xmin) *
            (canvas_width - 2L * plot_margin - 1L))
  py <- canvas_height - plot_margin -
    round((y - ymin) / (ymax - ymin) *
            (canvas_height - 2L * plot_margin - 1L))
  cbind(px = as.integer(px), py = as.integer(py))
}

set_pixel <- function(img, x, y, col) {
  ok <- x >= 1L & x <= ncol(img) & y >= 1L & y <= nrow(img)
  if (any(ok))
    img[cbind(y[ok], x[ok])] <- col
  img
}

# Integer line rasteriser.
draw_line <- function(img, x0, y0, x1, y1, col, width = 1L) {
  if (any(!is.finite(c(x0, y0, x1, y1))))
    return(img)

  n <- max(abs(x1 - x0), abs(y1 - y0)) + 1L
  if (n < 2L) n <- 2L
  xs <- round(seq(x0, x1, length.out = n))
  ys <- round(seq(y0, y1, length.out = n))

  offsets <- if (width <= 1L) 0L else seq.int(-floor(width / 2), floor(width / 2))
  for (dx in offsets) {
    for (dy in offsets) {
      img <- set_pixel(img, xs + dx, ys + dy, col)
    }
  }
  img
}

# Filled circular point marker.
draw_disc <- function(img, x, y, radius, col, outline = NULL) {
  if (!is.finite(x) || !is.finite(y))
    return(img)

  r <- as.integer(radius)
  gx <- seq.int(-r, r)
  gy <- seq.int(-r, r)
  dd <- expand.grid(dx = gx, dy = gy)
  keep <- dd$dx^2 + dd$dy^2 <= r^2
  dd <- dd[keep, , drop = FALSE]
  img <- set_pixel(img, x + dd$dx, y + dd$dy, col)

  if (!is.null(outline) && r >= 2L) {
    rr <- dd$dx^2 + dd$dy^2
    edge <- rr >= (r - 1)^2
    img <- set_pixel(img, x + dd$dx[edge], y + dd$dy[edge], outline)
  }
  img
}

new_canvas <- function() {
  img <- matrix(COL_BACKGROUND,
                nrow = canvas_height,
                ncol = canvas_width)

  # 500-unit reference grid.
  for (v in seq(0, 3000, by = 500)) {
    a <- xy_to_pixel(v, ymin)
    b <- xy_to_pixel(v, ymax)
    img <- draw_line(img, a[1, 'px'], a[1, 'py'],
                     b[1, 'px'], b[1, 'py'], COL_GRID)

    a <- xy_to_pixel(xmin, v)
    b <- xy_to_pixel(xmax, v)
    img <- draw_line(img, a[1, 'px'], a[1, 'py'],
                     b[1, 'px'], b[1, 'py'], COL_GRID)
  }

  # Landscape border.
  ll <- xy_to_pixel(xmin, ymin)
  lr <- xy_to_pixel(xmax, ymin)
  ul <- xy_to_pixel(xmin, ymax)
  ur <- xy_to_pixel(xmax, ymax)
  img <- draw_line(img, ll[1,1], ll[1,2], lr[1,1], lr[1,2], COL_BORDER, 2)
  img <- draw_line(img, lr[1,1], lr[1,2], ur[1,1], ur[1,2], COL_BORDER, 2)
  img <- draw_line(img, ur[1,1], ur[1,2], ul[1,1], ul[1,2], COL_BORDER, 2)
  img <- draw_line(img, ul[1,1], ul[1,2], ll[1,1], ll[1,2], COL_BORDER, 2)
  img
}

###############################################################################
### 5. Build one raster frame from PointHistory and EventLog
###############################################################################

build_frame <- function(timestep) {
  img <- new_canvas()

  if (timestep == 0L) {
    pp <- data.frame(
      id = seq_len(nrow(initial_points)),
      x = initial_points$x,
      y = initial_points$y,
      stage = initial_points$stage
    )
    new_ids <- integer(0)
    ev <- fit$EventLog[0, , drop = FALSE]
  } else {
    pp <- fit$PointHistory[
      fit$PointHistory$timestep == timestep,
      c('id', 'x', 'y', 'stage'),
      drop = FALSE
    ]
    ev <- fit$EventLog[fit$EventLog$timestep == timestep, , drop = FALSE]
    new_ids <- ev$point_id[ev$event == 'establishment']
  }

  # Reproductive dispersal lines for this timestep.
  if (timestep > 0L) {
    repro <- ev[ev$event == 'establishment', , drop = FALSE]
    if (nrow(repro) > 0L) {
      for (i in seq_len(nrow(repro))) {
        src <- parent_location_for_reproduction(repro$parent_id[i], timestep)
        dst <- c(repro$x[i], repro$y[i])
        if (all(is.finite(c(src, dst)))) {
          p0 <- xy_to_pixel(src[1], src[2])
          p1 <- xy_to_pixel(dst[1], dst[2])
          img <- draw_line(img,
                           p0[1,1], p0[1,2],
                           p1[1,1], p1[1,2],
                           COL_REPRO_LINE, 1L)
        }
      }
    }

    # Stage-transition movement lines for this timestep.
    trans <- ev[
      ev$event == 'stage_transition' &
        ev$detail == 'transition_kernel',
      , drop = FALSE
    ]
    if (nrow(trans) > 0L) {
      for (i in seq_len(nrow(trans))) {
        src <- point_location_before(trans$point_id[i], timestep)
        dst <- c(trans$x[i], trans$y[i])
        if (all(is.finite(c(src, dst)))) {
          p0 <- xy_to_pixel(src[1], src[2])
          p1 <- xy_to_pixel(dst[1], dst[2])
          img <- draw_line(img,
                           p0[1,1], p0[1,2],
                           p1[1,1], p1[1,2],
                           COL_TRANSITION_LINE, 2L)
        }
      }
    }
  }

  # Draw established points. Newly established stage-1 points are shown as
  # propagules for their birth timestep, then become ordinary stage-1 points.
  if (nrow(pp) > 0L) {
    pix <- xy_to_pixel(pp$x, pp$y)

    ordinary_stage1 <- pp$stage == 1L & !(pp$id %in% new_ids)
    stage2 <- pp$stage == 2L
    propagules <- pp$id %in% new_ids

    if (any(ordinary_stage1)) {
      for (i in which(ordinary_stage1))
        img <- draw_disc(img, pix[i,1], pix[i,2], 4L,
                         COL_STAGE1, outline = COL_BORDER)
    }

    if (any(stage2)) {
      for (i in which(stage2))
        img <- draw_disc(img, pix[i,1], pix[i,2], 4L,
                         COL_STAGE2, outline = COL_BORDER)
    }

    if (any(propagules)) {
      for (i in which(propagules))
        img <- draw_disc(img, pix[i,1], pix[i,2], 5L,
                         COL_PROPAGULE, outline = COL_BORDER)
    }
  }

  # Small, purely graphical legend in the upper-left margin:
  # blue = stage 1, teal = stage 2, red = new propagule.
  legend_y <- 12L
  img <- draw_disc(img, 20L, legend_y, 4L, COL_STAGE1, COL_BORDER)
  img <- draw_disc(img, 38L, legend_y, 4L, COL_STAGE2, COL_BORDER)
  img <- draw_disc(img, 56L, legend_y, 5L, COL_PROPAGULE, COL_BORDER)

  # Graphical timestep progress bar across the upper margin. This avoids a
  # font dependency while still making frame progression visible.
  progress_x0 <- 85L
  progress_x1 <- canvas_width - 20L
  progress_y <- 12L
  img <- draw_line(img, progress_x0, progress_y,
                   progress_x1, progress_y, COL_GRID, 3L)
  if (timestep > 0L) {
    done_x <- round(progress_x0 +
                      timestep / 18 * (progress_x1 - progress_x0))
    img <- draw_line(img, progress_x0, progress_y,
                     done_x, progress_y, COL_BORDER, 3L)
  }

  img
}

###############################################################################
### 6. Minimal GIF89a writer implemented in base R
###
### The encoder deliberately emits literal pixel codes and resets the LZW
### dictionary every 254 pixels. This keeps the GIF code width fixed at 9 bits,
### making the writer compact, transparent, and dependency-free.
###############################################################################

u16le <- function(x) {
  x <- as.integer(x)
  as.raw(c(x %% 256L, (x %/% 256L) %% 256L))
}

pack_9bit_codes <- function(codes) {
  codes <- as.integer(codes)
  n <- length(codes)
  nbytes <- as.integer(ceiling(n * 9 / 8))
  out <- integer(nbytes)

  pos <- (seq_len(n) - 1L) * 9L
  byte <- pos %/% 8L + 1L
  shift <- pos %% 8L

  # 9-bit values span no more than two output bytes.
  for (s in 0:7) {
    k <- which(shift == s)
    if (length(k) == 0L) next
    v <- bitwShiftL(codes[k], s)
    b <- byte[k]
    out[b] <- bitwOr(out[b], bitwAnd(v, 255L))
    hi <- bitwShiftR(v, 8L)
    b2 <- b + 1L
    ok <- b2 <= nbytes & hi != 0L
    if (any(ok))
      out[b2[ok]] <- bitwOr(out[b2[ok]], hi[ok])
  }

  as.raw(out)
}

gif_lzw_literal <- function(pixels) {
  pixels <- as.integer(pixels)
  if (any(pixels < 0L | pixels > 255L))
    stop('GIF pixels must be palette indices 0..255.')

  clear_code <- 256L
  end_code <- 257L
  chunk_size <- 254L
  n <- length(pixels)
  nchunks <- ceiling(n / chunk_size)

  codes <- integer(n + nchunks + 1L)
  pos <- 1L
  start <- 1L

  while (start <= n) {
    finish <- min(n, start + chunk_size - 1L)
    codes[pos] <- clear_code
    pos <- pos + 1L
    len <- finish - start + 1L
    codes[pos:(pos + len - 1L)] <- pixels[start:finish]
    pos <- pos + len
    start <- finish + 1L
  }

  codes[pos] <- end_code
  pack_9bit_codes(codes[seq_len(pos)])
}

write_subblocks <- function(con, raw_data) {
  n <- length(raw_data)
  start <- 1L
  while (start <= n) {
    len <- min(255L, n - start + 1L)
    writeBin(as.raw(len), con)
    writeBin(raw_data[start:(start + len - 1L)], con)
    start <- start + len
  }
  writeBin(as.raw(0L), con)
}

write_native_gif <- function(frames,
                             filename,
                             palette,
                             delay_cs = 100L,
                             loop = 0L) {
  if (length(frames) == 0L)
    stop('No frames supplied.')

  h <- nrow(frames[[1]])
  w <- ncol(frames[[1]])
  if (any(vapply(frames, nrow, integer(1)) != h) ||
      any(vapply(frames, ncol, integer(1)) != w))
    stop('All GIF frames must have the same dimensions.')

  con <- file(filename, open = 'wb')
  on.exit(close(con), add = TRUE)

  # Header and Logical Screen Descriptor.
  writeBin(charToRaw('GIF89a'), con)
  writeBin(u16le(w), con)
  writeBin(u16le(h), con)
  writeBin(as.raw(0xF7), con) # global table, 8-bit colour resolution, 256 entries
  writeBin(as.raw(COL_BACKGROUND), con)
  writeBin(as.raw(0L), con)

  # Global Colour Table: 256 RGB entries.
  writeBin(as.raw(as.vector(t(palette))), con)

  # NETSCAPE2.0 loop extension; loop = 0 means loop indefinitely.
  writeBin(as.raw(c(0x21, 0xFF, 0x0B)), con)
  writeBin(charToRaw('NETSCAPE2.0'), con)
  writeBin(as.raw(c(0x03, 0x01)), con)
  writeBin(u16le(loop), con)
  writeBin(as.raw(0x00), con)

  for (i in seq_along(frames)) {
    # Graphic Control Extension: disposal method 1, no transparency.
    writeBin(as.raw(c(0x21, 0xF9, 0x04, 0x04)), con)
    writeBin(u16le(delay_cs), con)
    writeBin(as.raw(c(0x00, 0x00)), con)

    # Image Descriptor.
    writeBin(as.raw(0x2C), con)
    writeBin(u16le(0L), con) # left
    writeBin(u16le(0L), con) # top
    writeBin(u16le(w), con)
    writeBin(u16le(h), con)
    writeBin(as.raw(0x00), con) # use global colour table

    # LZW minimum code size 8 (256-entry palette).
    writeBin(as.raw(0x08), con)

    # Flatten matrix in scanline order: left->right, top->bottom.
    pix <- as.integer(as.vector(t(frames[[i]])))
    compressed <- gif_lzw_literal(pix)
    write_subblocks(con, compressed)
  }

  writeBin(as.raw(0x3B), con) # GIF trailer
  invisible(filename)
}

###############################################################################
### 7. Build frames and write 1-second-per-timestep GIF
###############################################################################

frame_timesteps <- 0:18
frames <- lapply(frame_timesteps, build_frame)

output_gif <- file.path(demo_dir, 'illustrative_unmanaged_2stage_nativeR.gif')
write_native_gif(
  frames = frames,
  filename = output_gif,
  palette = palette_rgb,
  delay_cs = 100L,
  loop = 0L
)

# Animation-count table for interpretation/checking.
count_tab <- do.call(rbind, lapply(1:18, function(tt) {
  pp <- fit$PointHistory[fit$PointHistory$timestep == tt, , drop = FALSE]
  new_ids <- fit$EventLog$point_id[
    fit$EventLog$timestep == tt & fit$EventLog$event == 'establishment'
  ]
  data.frame(
    timestep = tt,
    total = nrow(pp),
    stage1 = sum(pp$stage == 1L & !(pp$id %in% new_ids)),
    stage2 = sum(pp$stage == 2L),
    propagules = sum(pp$id %in% new_ids),
    reproductive_dispersal_lines = sum(
      fit$EventLog$timestep == tt & fit$EventLog$event == 'establishment'
    ),
    transition_movement_lines = sum(
      fit$EventLog$timestep == tt &
        fit$EventLog$event == 'stage_transition' &
        fit$EventLog$detail == 'transition_kernel'
    )
  )
}))
write.csv(count_tab,
          file.path(output_dir, 'illustrative_unmanaged_timestep_counts_nativeR.csv'),
          row.names = FALSE)

cat('Native-R unmanaged 2-stage GIF complete.\n')
cat('GIF:', output_gif, '\n')
cat('Frames:', length(frames), '(timestep 0 plus 18 model timesteps)\n')
cat('Frame delay:', 1, 'second\n')
cat('Final population:', nrow(fit$FinalPoints), '\n')
print(count_tab)

q('no')
