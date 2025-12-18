-- `:help lsp`

vim.lsp.enable({
    "jetls",
})

-- Diagnostics configuration
vim.diagnostic.config({
    virtual_text = {
        prefix = "●",
        spacing = 2,
    },
    signs = {
        text = {
            [vim.diagnostic.severity.ERROR] = "E",
            [vim.diagnostic.severity.WARN] = "W",
            [vim.diagnostic.severity.INFO] = "I",
            [vim.diagnostic.severity.HINT] = "H",
        },
    },
    float = {
        border = "rounded",
        source = true,
    },
    severity_sort = true,
})

vim.api.nvim_create_autocmd("LspAttach", {
    callback = function(ev)
        local function map(mode, lhs, rhs, desc)
            vim.keymap.set(mode, lhs, rhs, { desc = desc, buffer = ev.buf })
        end

        vim.opt_local.omnifunc = "v:lua.vim.lsp:omnifunc"

        -- Navigation
        map("n", "gd", vim.lsp.buf.definition, "Go to definition")
        map("n", "gD", vim.lsp.buf.declaration, "Go to declaration")
        map("n", "gr", vim.lsp.buf.references, "Go to references")
        map("n", "gi", vim.lsp.buf.implementation, "Go to implementation")
        map("n", "gT", vim.lsp.buf.type_definition, "Go to type definition")

        -- Information
        map("n", "K", vim.lsp.buf.hover, "Hover documentation")
        map("i", "<C-k>", vim.lsp.buf.signature_help, "Signature help")

        -- Actions
        map("n", "<leader>rn", vim.lsp.buf.rename, "Rename symbol")
        map({ "n", "v" }, "<leader>ca", vim.lsp.buf.code_action, "Code action")
        map("n", "<leader>f", function()
            vim.lsp.buf.format({ async = true })
        end, "Format buffer")

        -- Diagnostics
        map("n", "[d", function()
            vim.diagnostic.jump({ count = -1, float = true })
        end, "Previous diagnostic")
        map("n", "]d", function()
            vim.diagnostic.jump({ count = 1, float = true })
        end, "Next diagnostic")
        map("n", "<leader>e", vim.diagnostic.open_float, "Show diagnostic")
        map("n", "<leader>q", vim.diagnostic.setloclist, "Diagnostics to loclist")
    end,
})
