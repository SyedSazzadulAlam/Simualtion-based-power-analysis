# PV Diagnostic: Simple Simulation-Based Power Analysis
# -----------------------------------------------------
# Earlier exploratory power analysis based on participant-level
# Pre-to-Post change scores.
#
# This is a simpler planning model than the trial-level GLMM simulation.
# It simulates binomial accuracy scores for each participant and compares
# the Pre-to-Post changes between Bimodal and Unimodal groups using a
# two-sample t-test.
#
# This script preserves the assumptions and scenarios from the original
# Power Analysis.R file.


# Scenario 1: Pilot-based probabilities ------------------------------------

set.seed(123)

# Current pilot probabilities
p_B_pre  <- 0.625
p_B_post <- 0.729

p_U_pre  <- 0.542
p_U_post <- 0.458

trials <- 16
nsim <- 5000


power_sim <- function(n_per_group) {

  p_values <- replicate(
    nsim,
    {

      # Simulate participant accuracy over 16 trials
      B_pre <- rbinom(
        n_per_group,
        trials,
        p_B_pre
      ) / trials

      B_post <- rbinom(
        n_per_group,
        trials,
        p_B_post
      ) / trials

      U_pre <- rbinom(
        n_per_group,
        trials,
        p_U_pre
      ) / trials

      U_post <- rbinom(
        n_per_group,
        trials,
        p_U_post
      ) / trials


      # Individual Pre-to-Post changes
      change_B <- B_post - B_pre
      change_U <- U_post - U_pre


      # Compare change scores between groups
      t.test(
        change_B,
        change_U
      )$p.value
    }
  )

  mean(p_values < 0.05)
}


sample_sizes <- c(
  5, 8, 10, 12,
  15, 20, 25, 30
)

pilot_results <- data.frame(
  n_per_group = sample_sizes,
  power = sapply(
    sample_sizes,
    power_sim
  )
)

print(pilot_results)


# Scenario 2: More conservative assumed effects ----------------------------

set.seed(123)

p_B_pre  <- 0.60
p_B_post <- 0.66

p_U_pre  <- 0.60
p_U_post <- 0.56

trials <- 16
nsim <- 5000


power_sim_conservative <- function(n_per_group) {

  p_values <- replicate(
    nsim,
    {

      B_pre <- rbinom(
        n_per_group,
        trials,
        p_B_pre
      ) / trials

      B_post <- rbinom(
        n_per_group,
        trials,
        p_B_post
      ) / trials

      U_pre <- rbinom(
        n_per_group,
        trials,
        p_U_pre
      ) / trials

      U_post <- rbinom(
        n_per_group,
        trials,
        p_U_post
      ) / trials


      change_B <- B_post - B_pre
      change_U <- U_post - U_pre


      t.test(
        change_B,
        change_U
      )$p.value
    }
  )

  mean(p_values < 0.05)
}


conservative_results <- data.frame(
  n_per_group = sample_sizes,
  power = sapply(
    sample_sizes,
    power_sim_conservative
  )
)

print(conservative_results)


# Optional paired-t-test reference calculation -----------------------------

# Original reference calculation:
# effect size d = 0.8
# alpha = .05
# target power = .80

if (requireNamespace("pwr", quietly = TRUE)) {

  paired_t_power <- pwr::pwr.t.test(
    d = 0.8,
    sig.level = 0.05,
    power = 0.80,
    type = "paired"
  )

  print(paired_t_power)

} else {

  message(
    "Optional package 'pwr' is not installed; ",
    "skipping the paired-t-test reference calculation."
  )
}
