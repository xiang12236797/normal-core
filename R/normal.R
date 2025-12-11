## normal.R - R wrappers for C-implemented Normal distribution

# 动态加载 C 编译出的 DLL
load_normal_dll <- function(path = NULL) {
  # 默认：在当前工作目录的 src 子目录下找 normal.dll / .so
  if (is.null(path)) {
    path <- file.path(getwd(), "src", paste0("normal", .Platform$dynlib.ext))
  }

  if (!file.exists(path)) {
    stop("Cannot find dynamic library at: ", path)
  }

  # 如果已经加载过，就不用重复加载
  if (!is.loaded("normal_pdf_c")) {
    dyn.load(path)
  }

  invisible(TRUE)
}

# ------------------------------
#   C 实现的 pdf / cdf / prob
# ------------------------------

# PDF 封装（类似于 stats::dnorm）
dnorm_c <- function(x, mean = 0, sd = 1) {
  x <- as.double(x)
  n <- as.integer(length(x))

  res <- .C("normal_pdf_c",
            n     = n,
            x     = x,
            mu    = as.double(mean),
            sigma = as.double(sd),
            out   = double(n))

  res$out
}

# CDF 封装（类似于 stats::pnorm）
pnorm_c <- function(q, mean = 0, sd = 1) {
  q <- as.double(q)
  n <- as.integer(length(q))

  res <- .C("normal_cdf_c",
            n     = n,
            x     = q,
            mu    = as.double(mean),
            sigma = as.double(sd),
            out   = double(n))

  res$out
}

# 区间概率：P(lower <= X <= upper)
pnorm_interval_c <- function(lower, upper,
                             mean = 0, sd = 1) {
  lower <- as.double(lower)
  upper <- as.double(upper)

  if (length(lower) != length(upper)) {
    stop("lower 和 upper 的长度必须相同")
  }

  n <- as.integer(length(lower))

  res <- .C("normal_prob_c",
            n     = n,
            lower = lower,
            upper = upper,
            mu    = as.double(mean),
            sigma = as.double(sd),
            out   = double(n))

  res$out
}
