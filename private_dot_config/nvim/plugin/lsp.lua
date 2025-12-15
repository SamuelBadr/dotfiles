-- `:help lsp`

vim.lsp.enable({
    "jetls",
})

vim.api.nvim_create_autocmd("LspAttach", {
    callback = function()
        -- Keymaps
        local function nmap(lhs, rhs, desc)
            vim.keymap.set("n", lhs, rhs, { desc = desc, buffer = 0 })
        end
        vim.opt_local.omnifunc = "v:lua.vim.lsp:omnifunc"
        nmap("gT", vim.lsp.buf.type_definition, "LSP Type definition")
    end,
})
