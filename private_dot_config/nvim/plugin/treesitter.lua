vim.pack.add({
    {
        src = "https://github.com/nvim-treesitter/nvim-treesitter",
        version = "main",
    },
})

vim.api.nvim_create_autocmd("PackChanged", {
    group = vim.api.nvim_create_augroup("myconfig.treesitter", { clear = true }),
    pattern = { "nvim-treesitter" },
    callback = function()
        vim.notify("Updating treesitter parsers", vim.log.levels.INFO)
        require("nvim-treesitter").update({ summary = true }):wait(30 * 1000)
    end,
})

require("nvim-treesitter").install({
    "julia",
    "zsh",
    "lua",
})

vim.api.nvim_create_autocmd("FileType", {
    pattern = require("nvim-treesitter").get_installed(),
    callback = function(ev)
        pcall(vim.treesitter.start, ev.buf)
    end,
})

local links = {
  ["@function"]        = "Function",
  ["@method"]          = "Function",
  ["@keyword"]         = "Keyword",
  ["@variable"]        = "Identifier",
  ["@parameter"]       = "Identifier",
  ["@field"]           = "Identifier",
  ["@property"]        = "Identifier",
  ["@string"]          = "String",
  ["@string.special"]  = "Special",
  ["@number"]          = "Number",
  ["@boolean"]         = "Boolean",
  ["@type"]            = "Type",
  ["@constant"]        = "Constant",
  ["@comment"]         = "Comment",
  ["@operator"]        = "Operator",
  ["@punctuation"]     = "Delimiter",
}

for capture, group in pairs(links) do
  vim.api.nvim_set_hl(0, capture, { link = group })
end

