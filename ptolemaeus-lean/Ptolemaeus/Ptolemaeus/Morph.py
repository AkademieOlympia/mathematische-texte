import asyncio
from morph_sandbox import MorphSandbox

async def main():
    
   # Use context manager for automatic cleanup
   async with await MorphSandbox.create() as sandbox:

      # Execute Python code directly
      result = await sandbox.execute_code("x = 42")

      result = await sandbox.execute_code("print(f'the answer is {x}')")
      print(result["output"])

asyncio.run(main())