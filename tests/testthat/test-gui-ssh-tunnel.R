test_that("SSH tunnel settings are validated and translated without a shell", {
  identity <- tempfile("liber-ssh-key-")
  writeLines("test fixture", identity)
  config <- list(
    host = "login.cluster.example.org", user = "researcher", port = 22L,
    remote_host = "127.0.0.1", remote_port = 8000L, local_port = 49199L,
    proxy_host = "gateway.example.org", proxy_user = "researcher",
    proxy_port = 22L, identity_file = identity,
    accept_new_host_key = TRUE, auto_start = TRUE
  )
  command <- LibeRation:::.liber_ssh_tunnel_command(config, ssh = "ssh")
  on.exit(unlink(command$ssh_config_file, force = TRUE), add = TRUE)
  expect_equal(command$local_port, 49199L)
  expect_equal(command$url, "http://127.0.0.1:49199")
  expect_true("127.0.0.1:49199:127.0.0.1:8000" %in% command$args)
  expect_true("StrictHostKeyChecking=accept-new" %in% command$args)
  expect_true("-F" %in% command$args)
  expect_identical(tail(command$args, 1L), "liber-managed-target")
  ssh_config <- readLines(command$ssh_config_file, warn = FALSE)
  expect_true(any(grepl("HostName gateway.example.org", ssh_config, fixed = TRUE)))
  expect_true(any(grepl("User researcher", ssh_config, fixed = TRUE)))
  expect_true(any(grepl(
    normalizePath(identity, winslash = "/"), ssh_config, fixed = TRUE
  )))
  expect_true("  ProxyJump liber-managed-gateway" %in% ssh_config)
  expect_false(any(grepl("StrictHostKeyChecking no", ssh_config, fixed = TRUE)))

  agent_command <- LibeRation:::.liber_ssh_tunnel_command(
    utils::modifyList(config, list(identity_file = "")), ssh = "ssh"
  )
  expect_true("researcher@gateway.example.org" %in% agent_command$args)
  expect_true("-J" %in% agent_command$args)
  expect_identical(agent_command$ssh_config_file, "")

  expect_error(
    LibeRation:::.liber_ssh_tunnel_normalize(
      utils::modifyList(config, list(host = "-oProxyCommand=bad"))
    ),
    "SSH host"
  )
  expect_error(
    LibeRation:::.liber_ssh_tunnel_normalize(
      utils::modifyList(config, list(local_port = 70000L))
    ),
    "local forwarded port"
  )
})

test_that("SSH remote definitions survive client-settings upgrades", {
  workspace <- nm_workspace(tempfile("liber-ssh-settings-"))
  remote <- list(
    name = "Research cluster", connection_mode = "ssh_tunnel", token = "secret",
    timeout = 30, user = "researcher", url = "http://127.0.0.1:49199",
    ssh = list(
      host = "login.cluster.example.org", user = "researcher", port = 22L,
      remote_host = "127.0.0.1", remote_port = 8000L, local_port = 0L,
      proxy_host = "", proxy_user = "", proxy_port = 22L,
      identity_file = "", accept_new_host_key = TRUE, auto_start = TRUE
    )
  )
  LibeRation:::.liber_client_settings_write(
    workspace, selected_queue = "cluster", remotes = list(cluster = remote)
  )
  restored <- LibeRation:::.liber_client_settings_read(workspace)
  expect_equal(restored$version, 9L)
  expect_equal(restored$selected_queue, "cluster")
  expect_equal(restored$remotes$cluster$connection_mode, "ssh_tunnel")
  expect_equal(restored$remotes$cluster$ssh$host, "login.cluster.example.org")
  expect_equal(restored$remotes$cluster$ssh$local_port, 0L)
})

test_that("SSH readiness discovers local keys without reading private material", {
  home <- tempfile("liber-ssh-home-")
  directory <- file.path(home, ".ssh")
  dir.create(directory, recursive = TRUE)
  private <- file.path(directory, "id_ed25519_liber")
  public_value <- paste(
    "ssh-ed25519",
    "AAAAC3NzaC1lZDI1NTE5AAAAIKbELKIjy0nYVnBIYqPCDZqOv6tsYSabGADBRuaMvMtF",
    "LibeR-test"
  )
  writeLines("private fixture is never parsed", private)
  writeLines(public_value, paste0(private, ".pub"))

  keys <- LibeRation:::.liber_ssh_key_candidates(home)
  expect_length(keys, 1L)
  expect_equal(keys[[1L]]$path, normalizePath(private, winslash = "/"))
  expect_true(keys[[1L]]$has_public_key)
  expect_equal(keys[[1L]]$public_key, public_value)
  expect_equal(LibeRation:::.liber_ssh_public_key(private), public_value)
})

test_that("SSH readiness reports the operating-system client when available", {
  ssh <- LibeRation:::.liber_ssh_tool("ssh")
  skip_if(!nzchar(ssh), "OpenSSH client is absent on this test runner")
  readiness <- LibeRation:::.liber_ssh_readiness()
  expect_true(readiness$available)
  expect_equal(readiness$ssh_path, ssh)
  expect_match(readiness$version, "OpenSSH", ignore.case = TRUE)
  expect_true(readiness$agent$status %in% c("ready", "empty", "unavailable"))

  installer <- paste(deparse(body(LibeRation:::.liber_ssh_install_client)),
                     collapse = "\n")
  expect_match(installer, "OpenSSH.Client~~~~0.0.1.0", fixed = TRUE)
  expect_match(installer, "Set-Service -Name 'ssh-agent'", fixed = TRUE)
  expect_match(installer, "Start-Service -Name 'ssh-agent'", fixed = TRUE)
})

test_that("the remote-server dialog exposes the managed SSH wizard", {
  source <- paste(readLines(
    system.file("htmlwidgets", "liberWorkbench.js", package = "LibeRation"),
    warn = FALSE
  ), collapse = "\n")
  expect_match(source, "SSH tunnel", fixed = TRUE)
  expect_match(source, "SSH destination", fixed = TRUE)
  expect_match(source, "Optional jump host", fixed = TRUE)
  expect_match(source, "queue_test", fixed = TRUE)
  expect_match(source, "Accept a new host key, but reject a changed key", fixed = TRUE)
  expect_match(source, "SSH readiness", fixed = TRUE)
  expect_match(source, "Install & enable OpenSSH", fixed = TRUE)
  expect_match(source, "Generate protected key", fixed = TRUE)
  expect_match(source, "Install on destination", fixed = TRUE)
  expect_match(source, "Test destination through gateway", fixed = TRUE)
  expect_match(source, "login.cluster.example.org", fixed = TRUE)
  expect_match(source, "gateway.example.org", fixed = TRUE)
  institutional_destination <- paste(c("myriad", "rc", "ucl", "ac", "uk"), collapse = ".")
  institutional_gateway <- paste(c("ssh-gateway", "ucl", "ac", "uk"), collapse = ".")
  expect_false(grepl(institutional_destination, source, fixed = TRUE))
  expect_false(grepl(institutional_gateway, source, fixed = TRUE))

  app <- liber_gui(workspace = tempfile("liber-ssh-gui-"), queue = FALSE,
                   launch.browser = NULL)
  server_function <- app[["serverFuncSource"]]()
  body_text <- paste(deparse(body(server_function)), collapse = "\n")
  expect_match(body_text, "start_managed_tunnel", fixed = TRUE)
  expect_match(body_text, "Managed SSH tunnels are unavailable", fixed = TRUE)
  expect_match(body_text, "connection_test", fixed = TRUE)
  expect_match(body_text, "ssh_install_client", fixed = TRUE)
  expect_match(body_text, "ssh_generate_key", fixed = TRUE)
  expect_match(body_text, "ssh_install_public_key", fixed = TRUE)
  expect_match(body_text, "ssh_test_hop", fixed = TRUE)
  expect_false(grepl(institutional_gateway, body_text, fixed = TRUE))
})
