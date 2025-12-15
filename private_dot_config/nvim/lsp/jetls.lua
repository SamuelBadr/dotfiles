-- https://aviatesk.github.io/JETLS.jl/release/#Neovim

return {
    cmd = {
        "jetls",
        "--threads=auto",
        "--",
    },
    filetypes = {"julia"},
    root_markers = {
        "Project.toml",
        ".git",
    },
}

