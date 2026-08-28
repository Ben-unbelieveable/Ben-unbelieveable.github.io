"""进阶 MCP Server：Tools + Resources + Prompts +（可选）Tasks。"""
from __future__ import annotations

import asyncio
import json
from typing import Any

from mcp.server.fastmcp import Context, FastMCP

APP_CONFIG: dict[str, Any] = {
    "app_name": "mcp-advanced-demo",
    "version": "1.0.0",
    "max_batch_chunks": 10,
}

mcp = FastMCP("advanced_demo")


@mcp.tool()
def count_words(text: str) -> dict[str, int]:
    """统计文本词数（空白分词，简单演示）。

    Args:
        text: 输入文本

    Returns:
        含 word_count 的字典
    """
    words = [w for w in text.split() if w.strip()]
    return {"word_count": len(words)}


@mcp.resource("config://app")
def app_config() -> str:
    """应用只读配置（JSON）。"""
    return json.dumps(APP_CONFIG, ensure_ascii=False, indent=2)


@mcp.prompt()
def code_review(language: str, focus: str = "security") -> str:
    """生成代码审查提示模板。

    Args:
        language: 编程语言
        focus: 审查重点，如 security / performance
    """
    return (
        f"请对以下 {language} 代码进行审查，重点关注 {focus}。"
        "列出问题、风险等级与修改建议。"
    )


@mcp.tool(task=True)
async def batch_process(chunks: int, ctx: Context) -> str:
    """模拟分块批处理（后台 Task + 进度）。

    Args:
        chunks: 分块数量，1–10

    Returns:
        完成摘要
    """
    n = max(1, min(chunks, APP_CONFIG["max_batch_chunks"]))
    for i in range(n):
        await asyncio.sleep(0.3)
        await ctx.report_progress(i + 1, n, f"chunk {i + 1}/{n}")
    return f"processed {n} chunks"


if __name__ == "__main__":
    import sys

    transport = "stdio"
    if len(sys.argv) > 1 and sys.argv[1] == "--http":
        transport = "streamable-http"
    mcp.run(transport=transport)
