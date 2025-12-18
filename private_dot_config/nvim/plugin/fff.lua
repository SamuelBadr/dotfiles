vim.pack.add({ "https://github.com/dmtrKovalenko/fff.nvim" })

-- Automatically download or build binary when plugin is installed/updated
vim.api.nvim_create_autocmd("PackChanged", {
    group = vim.api.nvim_create_augroup("myconfig.fff", { clear = true }),
    pattern = { "fff.nvim" },
    callback = function(ev)
        if ev.data.kind == "install" or ev.data.kind == "update" then
            vim.schedule(function()
                require("fff.download").download_or_build_binary()
            end)
        end
    end,
})

-- Plugin configuration
vim.g.fff = {
    lazy_sync = true,
}

-- Keymaps
vim.keymap.set("n", "<leader>ff", function()
    require("fff").find_files()
end, { desc = "Find files" })

vim.keymap.set("n", "<leader>fg", function()
    require("fff").find_in_git_root()
end, { desc = "Find files in git root" })

vim.keymap.set("n", "gf", function()
    require("fff").open_file_under_cursor()
end, { desc = "Open file under cursor" })
