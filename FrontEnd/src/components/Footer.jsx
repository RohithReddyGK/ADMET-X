import { FaReact } from "react-icons/fa";
import { SiTailwindcss, SiFlask, SiPython } from "react-icons/si";
import { MdScience } from "react-icons/md";
import { TbDatabaseCog } from "react-icons/tb";
import { motion } from "framer-motion";

const Footer = () => {
  return (
    <motion.footer
      className="bg-blue-100 dark:bg-gray-800 text-gray-800 dark:text-white py-5 shadow-md"
      initial={{ opacity: 0, y: 50 }}
      whileInView={{ opacity: 1, y: 0 }}
      viewport={{ once: true }}
      transition={{ duration: 1 }}
    >
      <div className="max-w-6xl mx-auto px-6 flex flex-col sm:flex-row items-center justify-between gap-4">
        <p className="font-bold text-sm sm:text-base text-center sm:text-left dark:text-gray-300">
          © {new Date().getFullYear()} ADMET-X Platform — All rights reserved.
        </p>
        <div className="flex flex-wrap justify-center sm:justify-start gap-4 text-4xl sm:text-5xl">
          <FaReact className="text-cyan-500 dark:text-cyan-400 hover:scale-110 transition" title="React" />
          <SiTailwindcss className="text-sky-400 dark:text-sky-400 hover:scale-110 transition" title="Tailwind CSS" />
          <SiFlask className="text-black dark:text-white hover:scale-110 transition" title="Flask" />
          <MdScience className="text-green-600 dark:text-green-400 hover:scale-110 transition" title="RDKit" />
          <TbDatabaseCog className="text-purple-600 dark:text-purple-400 hover:scale-110 transition" title="TDC" />
        </div>
      </div>
    </motion.footer>
  );
};

export default Footer;
